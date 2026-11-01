"""dadi wrappers for two-locus spectra and for triallelic/two-locus spectrum files (pack P22).

The functions work on spectra that are already computed: transformations (mask, fold, unfold,
project, ancestral misidentification), marginal single-locus spectra, mean r^2, linkage
statistics per bin, Poisson and multinomial log-likelihoods, optimal scaling, conditioning on a
derived count and dadi's text file format. A two-locus spectrum is a cube
``spectrum[AB][Ab][aB]`` of haplotype counts for ``n`` sampled chromosomes (ab is ``n`` minus the
three), given as nested lists with ``None`` for masked entries.

dadi's two-locus *integration* (``demographics.equilibrium_phi``/``equilibrium``/``two_epoch``/
``three_epoch``, ``integration.advance`` and the three optimizers that call it) is not offered:
it raises ``AttributeError: module 'dadi' has no attribute 'tridiag'`` for every input, and even
when that name is supplied it returns NaN (the C tridiagonal solver receives non-contiguous array
views). The genotype-spectrum functions (``genotype_spectrum_from_F`` and relatives,
``project_Gdict``, ``misidentification_genotype_dict``, ``genotype_exp_data_to_arrays``) were
written for Python 2 and fail for every input (``n / 2`` as an array size, ``dict.keys()[0]``).

dadi and NumPy are imported only when a function runs. Importing ``dadi.TwoLocus`` makes dadi create
its cache directory ``~/.dadi/TwoLocus_cache``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.population.dadi_p22_common import (
    choice,
    dadi_and_numpy,
    existing_file,
    flag,
    guarded,
    integer,
    masked_array,
    nested,
    new_file,
    number,
    result,
    scalar,
)

_TL_TRANSFORMS = (
    "mask_infeasible",
    "unfold",
    "fold",
    "fold_ancestral",
    "fold_lr",
    "project",
    "misidentification",
    "array_to_spectrum",
)
_FILE_KINDS = ("triallelic", "two_locus")
_MARGINAL_METHODS = ("marginal", "single_locus")
_R2_METHODS = ("numerics", "spectrum_method")
_MAX_TL_SAMPLE = 40
_MAX_TL_GRID = 40
_MAX_FILE_PRECISION = 17


def _cube(np: Any, value: object, name: str) -> Any:
    array = masked_array(np, value, name, (3,), nonnegative=False, square=True)
    if array.shape[0] - 1 > _MAX_TL_SAMPLE:
        raise ConfigurationError(f"{name} has more than {_MAX_TL_SAMPLE} sampled chromosomes.")
    return array


def _tl_spectrum(np: Any, value: object, name: str, folded: bool, infeasible: bool = False) -> Any:
    from dadi.TwoLocus import TLSpectrum_mod

    array = _cube(np, value, name)
    return TLSpectrum_mod.TLSpectrum(
        np.ma.getdata(array),
        mask=np.ma.getmaskarray(array),
        mask_infeasible=infeasible,
        data_folded=folded,
    )


def _tl_description(np: Any, fs: Any) -> dict[str, Any]:
    return {
        "sample_size": int(fs.shape[0]) - 1,
        "folded": bool(getattr(fs, "folded", False)),
        "spectrum": nested(np, fs),
    }


def two_locus_spectrum_transform(
    spectrum: Sequence[Sequence[Sequence[float | None]]],
    operation: str,
    *,
    folded: bool = False,
    sample_size: int | None = None,
    p_misid: float | None = None,
    apply_infeasible_mask: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask, fold, unfold, project or misidentify a two-locus spectrum, using ``dadi.TwoLocus``.

    ``operation`` is ``mask_infeasible`` (``TLSpectrum.mask_infeasible``), ``unfold``
    (``TLSpectrum.unfold``, needs ``folded=True``), ``fold`` (``TLSpectrum.fold``: fold when the
    ancestral state is unknown), ``fold_ancestral`` (``numerics.fold_ancestral``: track minor
    frequencies), ``fold_lr`` (``numerics.fold_lr``: fold by left/right allele), ``project``
    (``numerics.project`` to ``sample_size`` chromosomes), ``misidentification``
    (``numerics.misidentification`` with probability ``p_misid``) or ``array_to_spectrum``
    (``numerics.array_to_spectrum``: the values of ``spectrum`` with dadi's infeasible-entry mask).
    ``spectrum`` is an ``(n+1)^3`` nested list with ``None`` for masked entries; ``folded`` states
    whether it is folded. The input is not modified.
    """

    _dadi, np = dadi_and_numpy()
    name = choice(operation, _TL_TRANSFORMS, "operation")
    from dadi.TwoLocus import numerics

    fs = _tl_spectrum(
        np,
        spectrum,
        "spectrum",
        flag(folded, "folded"),
        flag(apply_infeasible_mask, "apply_infeasible_mask"),
    )
    extra: dict[str, Any] = {}
    if name == "mask_infeasible":
        out = guarded("mask_infeasible", fs.mask_infeasible)
    elif name == "unfold":
        out = guarded("unfold", fs.unfold)
    elif name == "fold":
        out = guarded("fold", fs.fold)
    elif name == "fold_ancestral":
        out = guarded("fold_ancestral", lambda: numerics.fold_ancestral(fs))
    elif name == "fold_lr":
        out = guarded("fold_lr", lambda: numerics.fold_lr(fs.copy()))
    elif name == "project":
        if sample_size is None:
            raise ConfigurationError("sample_size is required to project.")
        target = integer(sample_size, "sample_size", minimum=3, maximum=int(fs.shape[0]) - 1)
        extra["projected_to"] = target
        out = guarded("project", lambda: numerics.project(fs, target))
    elif name == "misidentification":
        if p_misid is None:
            raise ConfigurationError("p_misid is required for misidentification.")
        probability = number(p_misid, "p_misid", minimum=0.0, maximum=1.0)
        extra["p_misid"] = probability
        out = guarded("misidentification", lambda: numerics.misidentification(fs, probability))
    else:
        out = guarded("array_to_spectrum", lambda: numerics.array_to_spectrum(np.ma.getdata(fs)))
    return result(
        "two_locus_spectrum_transform",
        "spectrum[i][j][k] is the count of AB, Ab and aB haplotypes among n chromosomes; masked "
        "entries are None",
        {"transform": name, **extra, **_tl_description(np, out)},
    )


def two_locus_marginal_spectra(
    spectrum: Sequence[Sequence[Sequence[float | None]]],
    *,
    method: str = "marginal",
    folded: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Single-locus spectra of locus A and locus B from a two-locus spectrum, using dadi.

    ``method`` is ``marginal`` (``TLSpectrum.marginalA`` and ``marginalB``, which sum all entries,
    masked or not) or ``single_locus`` (``numerics.to_single_locus``, which skips masked entries).
    Entry ``k`` of each result is the number of sites with ``k`` derived alleles at that locus.
    """

    _dadi, np = dadi_and_numpy()
    way = choice(method, _MARGINAL_METHODS, "method")
    fs = _tl_spectrum(np, spectrum, "spectrum", flag(folded, "folded"))
    from dadi.TwoLocus import numerics

    if way == "marginal":
        a, b = guarded("marginal spectra", lambda: (fs.marginalA(), fs.marginalB()))
    else:
        a, b = guarded("single-locus spectra", lambda: numerics.to_single_locus(fs))
    return result(
        "two_locus_marginal_spectra",
        "marginal[k] is the expected number of sites with k derived alleles at the locus among n "
        "chromosomes",
        {
            "method": way,
            "sample_size": int(fs.shape[0]) - 1,
            "locus_a": nested(np, a),
            "locus_b": nested(np, b),
        },
    )


def two_locus_mean_r2(
    spectrum: Sequence[Sequence[Sequence[float | None]]],
    *,
    method: str = "numerics",
    folded: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean r^2 (squared correlation between the loci) of a two-locus spectrum, using dadi.

    ``method`` is ``numerics`` (``dadi.TwoLocus.numerics.mean_r2``, a loop over unmasked entries)
    or ``spectrum_method`` (``TLSpectrum.mean_r2``, which uses the r^2 of every bin). Both weight
    r^2 by the spectrum entries.
    """

    _dadi, np = dadi_and_numpy()
    way = choice(method, _R2_METHODS, "method")
    fs = _tl_spectrum(np, spectrum, "spectrum", flag(folded, "folded"))
    from dadi.TwoLocus import numerics

    value = guarded("mean r2", (lambda: numerics.mean_r2(fs)) if way == "numerics" else fs.mean_r2)
    return result(
        "two_locus_mean_r2",
        "spectrum-weighted mean of r^2 = D^2 / (pA (1 - pA) pB (1 - pB))",
        {"method": way, "sample_size": int(fs.shape[0]) - 1, "mean_r2": scalar(value)},
    )


def two_locus_ld_per_bin(sample_size: int, *, verbose: bool = True) -> dict[str, Any]:
    """Linkage disequilibrium D and r^2 for every haplotype-count bin.

    Uses ``numerics.LD_per_bin``.

    Entry ``[AB][Ab][aB]`` of ``D`` is ``pAB - pA pB`` and of ``r2`` is ``D^2 / (pA (1 - pA) pB
    (1 - pB))`` for ``sample_size`` chromosomes; undefined bins are ``None``.
    """

    _dadi, np = dadi_and_numpy()
    n = integer(sample_size, "sample_size", minimum=2, maximum=_MAX_TL_SAMPLE)
    from dadi.TwoLocus import numerics

    d_bins, r2_bins = guarded("LD_per_bin", lambda: numerics.LD_per_bin(n))
    return result(
        "two_locus_ld_per_bin",
        "D = pAB - pA pB; r2 = D^2 / (pA (1 - pA) pB (1 - pB)) with p the haplotype counts over n",
        {"sample_size": n, "D": nested(np, d_bins), "r2": nested(np, r2_bins)},
    )


def two_locus_log_likelihood(
    model: Sequence[Sequence[Sequence[float | None]]],
    data: Sequence[Sequence[Sequence[float | None]]],
    *,
    multinom: bool = False,
    per_bin: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Poisson or multinomial log-likelihood of a two-locus data spectrum.

    Uses ``dadi.TwoLocus.inference``.

    ``multinom=False`` is ``ll`` (the model is used as given); ``multinom=True`` first rescales the
    model by its optimal scalar (``ll_multinom``). ``per_bin=True`` returns the value of every bin
    (``ll_per_bin``, ``ll_multinom_per_bin``) instead of their sum. ``model`` and ``data`` are
    ``(n+1)^3`` nested lists of the same size, ``None`` for masked entries.
    """

    _dadi, np = dadi_and_numpy()
    multi, each = flag(multinom, "multinom"), flag(per_bin, "per_bin")
    expected = _tl_spectrum(np, model, "model", False)
    observed = _tl_spectrum(np, data, "data", False)
    if expected.shape != observed.shape:
        raise ConfigurationError("model and data must have the same shape.")
    from dadi.TwoLocus import inference

    if each:
        function = inference.ll_multinom_per_bin if multi else inference.ll_per_bin
    else:
        function = inference.ll_multinom if multi else inference.ll
    value = guarded("two-locus log-likelihood", lambda: function(expected, observed))
    data_out: dict[str, Any] = {
        "multinom": multi,
        "per_bin": each,
        "sample_size": int(expected.shape[0]) - 1,
    }
    if each:
        data_out["log_likelihood_per_bin"] = nested(np, value)
    else:
        data_out["log_likelihood"] = scalar(value)
    return result(
        "two_locus_log_likelihood",
        "Poisson log-likelihood of the data counts given the model",
        data_out,
    )


def two_locus_log_likelihood_over_rho_bins(
    models: Sequence[Sequence[Sequence[Sequence[float | None]]]],
    data: Sequence[Sequence[Sequence[Sequence[float | None]]]],
    *,
    multinom: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Total log-likelihood over recombination-rate bins, using ``inference.ll_over_rho_bins``.

    ``models`` and ``data`` are lists with one two-locus spectrum per rho bin. The result is the sum
    of the per-bin Poisson log-likelihoods, or with ``multinom=True`` of the multinomial ones
    (``ll_over_rho_bins_multinom``, each bin rescaled by its own optimal scalar).
    """

    _dadi, np = dadi_and_numpy()
    multi = flag(multinom, "multinom")
    for name, value in (("models", models), ("data", data)):
        if not isinstance(value, Sequence) or not value:
            raise ConfigurationError(f"{name} must be a non-empty list of two-locus spectra.")
    if len(models) != len(data):
        raise ConfigurationError("models and data must have one spectrum per rho bin each.")
    expected = [_tl_spectrum(np, item, "models entry", False) for item in models]
    observed = [_tl_spectrum(np, item, "data entry", False) for item in data]
    from dadi.TwoLocus import inference

    function = inference.ll_over_rho_bins_multinom if multi else inference.ll_over_rho_bins
    value = guarded("log-likelihood over rho bins", lambda: function(expected, observed))
    return result(
        "two_locus_log_likelihood_over_rho_bins",
        "sum over rho bins of the two-locus Poisson log-likelihood",
        {"multinom": multi, "num_bins": len(expected), "log_likelihood": scalar(value)},
    )


def scale_two_locus_spectrum(
    model: Sequence[Sequence[Sequence[float | None]]],
    data: Sequence[Sequence[Sequence[float | None]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Optimal scalar between a two-locus model and data, and the scaled model (dadi ``inference``).

    The scalar is the data total over the model total on the bins unmasked in both
    (``optimal_sfs_scaling``); the scaled model is ``scalar * model`` (``optimally_scaled_sfs``).
    """

    _dadi, np = dadi_and_numpy()
    expected = _tl_spectrum(np, model, "model", False)
    observed = _tl_spectrum(np, data, "data", False)
    if expected.shape != observed.shape:
        raise ConfigurationError("model and data must have the same shape.")
    from dadi.TwoLocus import inference

    scaling = guarded("optimal scaling", lambda: inference.optimal_sfs_scaling(expected, observed))
    scaled = guarded(
        "optimally scaled spectrum", lambda: inference.optimally_scaled_sfs(expected, observed)
    )
    return result(
        "scale_two_locus_spectrum",
        "scalar = sum(data) / sum(model) over the bins unmasked in both",
        {
            "scaling": scalar(scaling),
            "sample_size": int(expected.shape[0]) - 1,
            "scaled_model": nested(np, scaled),
        },
    )


def condition_two_locus_spectrum(
    spectrum: Sequence[Sequence[Sequence[float | None]]],
    derived_count: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Two-locus spectrum conditioned on ``derived_count`` derived alleles at locus A.

    dadi methods: ``numerics.condition``.

    Entry ``[i][j]`` of the result is the spectrum value with ``i`` AB haplotypes and ``j`` aB
    haplotypes, for sites where ``AB + Ab`` equals ``derived_count``. Its shape is
    ``(derived_count + 1) x (n - derived_count + 1)``.
    """

    _dadi, np = dadi_and_numpy()
    fs = _tl_spectrum(np, spectrum, "spectrum", False)
    n = int(fs.shape[0]) - 1
    count = integer(derived_count, "derived_count", minimum=1, maximum=n - 1)
    from dadi.TwoLocus import numerics

    conditioned = guarded("condition", lambda: numerics.condition(fs, count))
    return result(
        "condition_two_locus_spectrum",
        "conditioned[i][j] = spectrum[i][count - i][j]: i AB and j aB haplotypes, "
        "given count derived alleles at A",
        {"sample_size": n, "derived_count": count, "spectrum": nested(np, conditioned)},
    )


def sample_two_locus_density(
    phi: Sequence[Sequence[Sequence[float]]],
    sample_size: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected two-locus spectrum sampled from a density on the grid.

    Uses ``numerics.sample_cached``.

    ``phi`` is the density ``phi[AB][Ab][aB]`` on the uniform grid ``linspace(0, 1, g + 1)``
    (a cube of side ``g + 1``, at most 41); the result is the multinomial-weighted integral for
    ``sample_size`` chromosomes. dadi's own forward integration of the density is not offered (see
    the module notes).
    """

    _dadi, np = dadi_and_numpy()
    density = masked_array(np, phi, "phi", (3,), nonnegative=False, minimum_length=4, square=True)
    grid = int(density.shape[0]) - 1
    if grid > _MAX_TL_GRID:
        raise ConfigurationError(f"phi may have at most {_MAX_TL_GRID + 1} points per side.")
    n = integer(sample_size, "sample_size", minimum=2, maximum=_MAX_TL_SAMPLE)
    from dadi.TwoLocus import numerics

    x = numerics.grid(grid)
    spectrum = guarded(
        "sample_cached", lambda: numerics.sample_cached(np.ma.getdata(density), n, x)
    )
    return result(
        "sample_two_locus_density",
        "spectrum[i][j][k] is the expected count of AB, Ab and aB haplotypes among n chromosomes",
        {"grid_points": grid, **_tl_description(np, spectrum)},
    )


def _file_spectrum(
    np: Any,
    kind: str,
    value: object,
    folded_major: bool,
    folded_ancestral: bool,
    extrap_x: float | None,
    extrap_t: float | None,
) -> Any:
    if kind == "triallelic":
        from dadi.Triallele import TriSpectrum_mod

        array = masked_array(np, value, "spectrum", (2,), nonnegative=False, square=True)
        fs = TriSpectrum_mod.TriSpectrum(
            np.ma.getdata(array),
            mask=np.ma.getmaskarray(array),
            mask_infeasible=False,
            data_folded_major=folded_major,
            data_folded_ancestral=folded_ancestral,
            check_folding_major=False,
            check_folding_ancestral=False,
        )
    else:
        fs = _tl_spectrum(np, value, "spectrum", folded_major)
    fs.extrap_x = extrap_x
    fs.extrap_t = extrap_t
    return fs


def write_multilocus_spectrum_file(
    spectrum: Sequence[Any],
    destination: str,
    *,
    kind: str = "triallelic",
    folded_major: bool = False,
    folded_ancestral: bool = False,
    extrap_x: float | None = None,
    extrap_t: float | None = None,
    precision: int = 16,
    comment_lines: Sequence[str] = (),
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a triallelic or two-locus spectrum to a dadi text file.

    dadi methods: ``TriSpectrum.to_file``, ``TLSpectrum.to_file``.

    ``kind`` is ``triallelic`` (square ``(n+1)^2`` spectrum) or ``two_locus`` (``(n+1)^3`` cube;
    only ``folded_major`` is used as its folded flag). The file holds the comment lines, a header
    with the size, folding and extrapolation values, the data line and the mask line.
    An existing file is never overwritten.
    """

    _dadi, np = dadi_and_numpy()
    name = choice(kind, _FILE_KINDS, "kind")
    target = new_file(destination, "destination")
    digits = integer(precision, "precision", minimum=1, maximum=_MAX_FILE_PRECISION)
    comments = [str(line) for line in comment_lines]
    fs = _file_spectrum(
        np,
        name,
        spectrum,
        flag(folded_major, "folded_major"),
        flag(folded_ancestral, "folded_ancestral"),
        None if extrap_x is None else number(extrap_x, "extrap_x"),
        None if extrap_t is None else number(extrap_t, "extrap_t"),
    )
    guarded(
        "write spectrum file",
        lambda: fs.to_file(str(target), precision=digits, comment_lines=comments),
    )
    return result(
        "write_multilocus_spectrum_file",
        "dadi spectrum text format: header, data line, mask line (1 = masked)",
        {
            "kind": name,
            "destination": str(target),
            "bytes": target.stat().st_size,
            "sample_size": int(fs.shape[0]) - 1,
        },
    )


def read_multilocus_spectrum_file(
    path: str,
    *,
    kind: str = "triallelic",
    mask_infeasible: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a triallelic or two-locus spectrum from a dadi text file.

    dadi methods: ``TriSpectrum.from_file``, ``TLSpectrum.from_file``.

    ``kind`` must match the file. ``mask_infeasible`` masks the entries that cannot occur, as dadi
    does when it reads a file. The result lists the spectrum (``None`` for masked entries), the
    folding flags, the extrapolation values and the file's comment lines.
    """

    _dadi, np = dadi_and_numpy()
    name = choice(kind, _FILE_KINDS, "kind")
    source = existing_file(path, "path")
    infeasible = flag(mask_infeasible, "mask_infeasible")
    if name == "triallelic":
        from dadi.Triallele import TriSpectrum_mod as module

        fs, comments = guarded(
            "read spectrum file",
            lambda: module.TriSpectrum.from_file(
                str(source), mask_infeasible=infeasible, return_comments=True
            ),
        )
        folding = {
            "folded_major": bool(fs.folded_major),
            "folded_ancestral": bool(fs.folded_ancestral),
        }
    else:
        from dadi.TwoLocus import TLSpectrum_mod as module

        fs, comments = guarded(
            "read spectrum file",
            lambda: module.TLSpectrum.from_file(
                str(source), mask_infeasible=infeasible, return_comments=True
            ),
        )
        folding = {"folded": bool(fs.folded)}
    return result(
        "read_multilocus_spectrum_file",
        "dadi spectrum text format; masked entries are None",
        {
            "kind": name,
            "path": str(source),
            "sample_size": int(fs.shape[0]) - 1,
            **folding,
            "extrap_x": None if fs.extrap_x is None else float(fs.extrap_x),
            "extrap_t": None if fs.extrap_t is None else float(fs.extrap_t),
            "comments": list(comments),
            "spectrum": nested(np, fs),
        },
    )


__all__ = [
    "condition_two_locus_spectrum",
    "read_multilocus_spectrum_file",
    "sample_two_locus_density",
    "scale_two_locus_spectrum",
    "two_locus_ld_per_bin",
    "two_locus_log_likelihood",
    "two_locus_log_likelihood_over_rho_bins",
    "two_locus_marginal_spectra",
    "two_locus_mean_r2",
    "two_locus_spectrum_transform",
    "write_multilocus_spectrum_file",
]

_register_public_verbose(globals(), tuple(__all__))
