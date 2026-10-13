"""dadi wrappers for SNP and VCF data, frequency-spectrum operations and demographic models.

The functions here wrap dadi's data readers, spectrum readers and writers, ancestral
misidentification, bootstrap and pooling helpers, the demographic models for one, two and three
populations, the demographic-plus-selection models, the optimizer and residual helpers, and the
Fisher-information uncertainty. Every input is local. Model functions are selected by name from
fixed lists, so no arbitrary code is run. Results are spectra as nested lists, with corner
entries that dadi masks returned as None.

Output files are never overwritten, and dadi is imported only when a function runs.
"""

from __future__ import annotations

import importlib
import random
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population.backends import (
    _DADI_EXTRA,
    _finite,
    _integer,
    _number,
    _provenance,
    _text,
)

_MAX_1D_SAMPLE = 2000
_MAX_2D_SAMPLE = 200
_MAX_3D_SAMPLE = 40
_MAX_BOOTSTRAPS = 1000
_MAX_MAXITER = 10_000
_OPTIMIZERS = {
    "bfgs": "optimize_log",
    "lbfgsb": "optimize_log_lbfgsb",
    "nelder-mead": "optimize_log_fmin",
    "powell": "optimize_log_powell",
}
_RESIDUALS = {"anscombe": "Anscombe_Poisson_residual", "poisson": "linear_Poisson_residual"}
_TRI_ORDER_BASES = "CGTA"  # dadi's internal base order for tri_freq_dict_to_array
_MODEL_CONVENTION = (
    "spectrum[i] is the count for derived-allele frequency class i of n sampled chromosomes; "
    "entries dadi masks (for example the corners) are None"
)

# Model families. Each family has one representative (the first name) and the rest are variants
# selected by the model argument. Names are the attributes of the dadi module that is used.
_DEMOGRAPHICS_1D = (
    "two_epoch",
    "growth",
    "bottlegrowth",
    "bottlegrowth_1d",
    "three_epoch",
    "three_epoch_inbreeding",
)
_DEMOGRAPHICS_2D = (
    "split_mig",
    "no_mig",
    "sym_mig",
    "asym_mig",
    "anc_sym_mig",
    "anc_asym_mig",
    "sec_contact_sym_mig",
    "sec_contact_asym_mig",
    "no_mig_size",
    "sym_mig_size",
    "asym_mig_size",
    "anc_sym_mig_size",
    "anc_asym_mig_size",
    "sec_contact_sym_mig_size",
    "sec_contact_asym_mig_size",
    "sym_mig_twoepoch",
    "asym_mig_twoepoch",
    "sec_contact_sym_mig_three_epoch",
    "sec_contact_asym_mig_three_epoch",
    "sec_contact_sym_mig_size_three_epoch",
    "sec_contact_asym_mig_size_three_epoch",
    "vic_no_mig",
    "vic_anc_sym_mig",
    "vic_anc_asym_mig",
    "vic_sec_contact_sym_mig",
    "vic_sec_contact_asym_mig",
    "founder_nomig",
    "founder_sym",
    "founder_asym",
    "vic_no_mig_admix_early",
    "vic_no_mig_admix_late",
    "vic_two_epoch_admix",
    "founder_nomig_admix_early",
    "founder_nomig_admix_late",
    "founder_nomig_admix_two_epoch",
    "bottlegrowth_2d",
    "bottlegrowth",
    "bottlegrowth_split",
    "bottlegrowth_split_mig",
    "split_asym_mig",
    "split_delay_mig",
    "IM",
    "IM_pre",
)
_DEMOGRAPHICS_3D = (
    "split_nomig",
    "split_symmig_all",
    "split_symmig_adjacent",
    "refugia_adj_1",
    "refugia_adj_2",
    "refugia_adj_3",
    "ancmig_adj_3",
    "ancmig_adj_2",
    "ancmig_adj_1",
    "sim_split_no_mig",
    "sim_split_no_mig_size",
    "sim_split_sym_mig_all",
    "sim_split_sym_mig_adjacent",
    "sim_split_refugia_sym_mig_all",
    "sim_split_refugia_sym_mig_adjacent",
    "split_nomig_size",
    "ancmig_2_size",
    "sim_split_refugia_sym_mig_adjacent_size",
    "refugia_adj_2_var_sym",
    "refugia_adj_2_var_uni",
    "refugia_adj_3_var_sym",
    "refugia_adj_3_var_uni",
    "split_sym_mig_adjacent_var1",
    "split_uni_mig_adjacent_var1",
    "split_sym_mig_adjacent_var2",
    "split_uni_mig_adjacent_var2",
    "sim_split_sym_mig_adjacent_var",
    "sim_split_uni_mig_adjacent_var",
    "sim_split_refugia_sym_mig_adjacent_var",
    "sim_split_refugia_uni_mig_adjacent_var",
    "admix_origin_no_mig",
    "admix_origin_sym_mig_adj",
    "admix_origin_uni_mig_adj",
    "out_of_africa",
)
_DEMOGRAPHICS_SELECTION = (
    "split_mig_sel",
    "equil",
    "two_epoch_sel",
    "IM_pre_sel",
    "IM_pre_sel_single_gamma",
    "IM_sel",
    "IM_sel_single_gamma",
    "split_mig_sel_single_gamma",
    "split_asym_mig_sel",
    "split_asym_mig_sel_single_gamma",
    "split_delay_mig_sel",
    "split_delay_mig_sel_single_gamma",
    "three_epoch_sel",
    "bottlegrowth_2d_sel",
    "bottlegrowth_2d_sel_single_gamma",
    "bottlegrowth_split_sel",
    "bottlegrowth_split_sel_single_gamma",
    "bottlegrowth_split_mig_sel",
    "bottlegrowth_split_mig_sel_single_gamma",
    "growth_sel",
    "bottlegrowth_1d_sel",
)
_FIT_OPTIMIZER_METHODS = tuple(_OPTIMIZERS)


def _new_file(destination: object, name: str) -> Path:
    if not isinstance(destination, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    target = Path(_text(str(destination), name)).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    return target


def _choice(value: object, options: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {sorted(options)}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _sizes(value: object, dims: int, maximum: int, name: str = "sample_sizes") -> tuple[int, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != dims:
        raise ConfigurationError(f"{name} must list {dims} integer sample sizes.")
    return tuple(_integer(item, "sample size", minimum=2, maximum=maximum) for item in value)


def _numbers(value: object, name: str, *, minimum: int = 1) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) < minimum:
        raise ConfigurationError(f"{name} must be a list of at least {minimum} numbers.")
    return [_number(item, f"{name} entry") for item in value]


def _bounds(value: object, count: int, name: str) -> list[float] | None:
    if value is None:
        return None
    values = _numbers(value, name)
    if len(values) != count:
        raise ConfigurationError(f"{name} must have {count} entries, one per parameter.")
    return values


def _load_snp(dadi: Any, path: Path) -> dict[str, Any]:
    try:
        return dadi.Misc.make_data_dict(str(path))
    except (ValueError, KeyError, IndexError, OSError) as exc:
        raise ConfigurationError(f"dadi could not parse {path.name} as an SNP file: {exc}") from exc


def _populations(data: Mapping[str, Any]) -> list[str]:
    """Population IDs in header order, taken from the calls of the first SNP."""

    for info in data.values():
        return list(info["calls"])
    return []


def _pop_ids(value: object, available: Sequence[str]) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("pop_ids must be a non-empty list of population IDs.")
    ids = [_text(item, "population ID") for item in value]
    if len(set(ids)) != len(ids):
        raise ConfigurationError("pop_ids must not repeat a population.")
    missing = [pop for pop in ids if pop not in available]
    if missing:
        raise ConfigurationError(
            f"pop_ids {missing} are not in the file; it has {list(available)}."
        )
    return ids


def _counts_array(
    np: Any,
    value: object,
    name: str,
    dims: tuple[int, ...] = (1, 2, 3),
    *,
    nonnegative: bool = True,
) -> Any:
    """A finite spectrum with 1 to 3 axes, each of length at least 3 (n >= 2).

    Observed counts must be non-negative; model spectra and linear transforms may be signed.
    """

    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular list of numbers.") from exc
    if array.ndim not in dims or min(array.shape) < 3:
        raise ConfigurationError(f"{name} must have {dims} axes, each of length at least 3.")
    if not np.isfinite(array).all() or (nonnegative and (array < 0).any()):
        raise ConfigurationError(f"{name} must contain only finite, non-negative values.")
    return array


def _nested(np: Any, value: Any) -> Any:
    """Nested lists of floats, with masked entries as None."""

    data = np.ma.getdata(value).tolist()
    mask = np.ma.getmaskarray(value).tolist()
    return _apply_mask(data, mask)


def _apply_mask(data: Any, mask: Any) -> Any:
    if isinstance(data, list):
        return [_apply_mask(item, flag) for item, flag in zip(data, mask, strict=True)]
    number = None if mask else _finite(data)
    return None if number is None else number + 0.0  # + 0.0 turns -0.0 into 0.0


def _spectrum_sizes(spectrum: Any) -> list[int]:
    return [int(n) for n in spectrum.sample_sizes]


def _run_model(function: Any, params: list[float], ns: tuple[int, ...], points: int) -> Any:
    try:
        return function(tuple(params), ns, points)
    except (ValueError, TypeError, IndexError, KeyError) as exc:
        raise ConfigurationError(
            f"dadi model {function.__name__} rejected params={params} for sample sizes {ns}: {exc}"
        ) from exc


def _model_function(family: Any, name: str, options: Sequence[str]) -> Any:
    model = _choice(name, options, "model")
    function = getattr(family, model, None)
    if function is None:
        raise ConfigurationError(f"Installed dadi does not provide model {model!r}.")
    return function


def read_snp_data(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read a dadi SNP file and summarize it, using dadi's ``Misc.make_data_dict()``.

    The file has a header line ``Context OutgroupContext Allele1 <pops...> Allele2 <id>`` and one
    row per SNP, with the reference-allele and derived-allele counts for each population. A SNP
    is polarized when its outgroup allele is one of the two segregating alleles.
    """

    dadi = require("dadi", _DADI_EXTRA)
    source = input_path(path)
    data = _load_snp(dadi, source)
    biallelic = [info for info in data.values() if len(info["segregating"]) == 2]
    polarized = [
        info
        for info in biallelic
        if info.get("outgroup_allele", "-") != "-"
        and info["outgroup_allele"] in info["segregating"]
    ]
    return {
        "status": "completed",
        "operation": "read_snp_data",
        "source": str(source),
        "populations": _populations(data),
        "num_snps": len(data),
        "num_biallelic": len(biallelic),
        "num_polarized": len(polarized),
        "provenance": _provenance("dadi", "dadi SNP file; biallelic SNPs, outgroup-polarized"),
    }


def snp_configuration_counts(
    path: str | Path,
    *,
    pop_ids: Sequence[str],
    verbose: bool = True,
) -> dict[str, Any]:
    """Count SNP configurations for the given populations, using dadi's ``Misc.count_data_dict()``.

    Each row gives the successful calls and derived calls per population, whether the SNP is
    polarized, and the number of SNPs with that configuration. Rows are sorted by configuration.
    Non-biallelic SNPs are skipped, as dadi does.
    """

    dadi = require("dadi", _DADI_EXTRA)
    source = input_path(path)
    data = _load_snp(dadi, source)
    selected = _pop_ids(pop_ids, _populations(data))
    counts = dadi.Misc.count_data_dict(data, selected)
    rows = [
        {
            "successful_calls": [int(n) for n in key[0]],
            "derived_calls": [int(n) for n in key[1]],
            "polarized": bool(key[2]),
            "count": int(number),
        }
        for key, number in sorted(counts.items(), key=lambda item: item[0])
    ]
    return {
        "status": "completed",
        "operation": "snp_configuration_counts",
        "source": str(source),
        "pop_ids": selected,
        "num_configurations": len(rows),
        "num_snps_counted": sum(row["count"] for row in rows),
        "configurations": rows,
        "provenance": _provenance("dadi", "dadi count_data_dict; biallelic SNPs only"),
    }


def read_vcf_data(
    vcf_path: str | Path,
    popinfo_path: str | Path,
    *,
    parser: str = "builtin",
    filter: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a VCF file and a population table into a dadi data dictionary, and summarize it.

    ``popinfo_path`` has one ``sample population`` pair per line. ``parser="builtin"`` uses
    dadi's own VCF reader; ``parser="cyvcf2"`` uses cyvcf2, which needs the ``cyvcf2`` package.
    Both give the same dictionary for the same file. ``filter=True`` skips filtered sites.
    """

    dadi = require("dadi", _DADI_EXTRA)
    mode = _choice(parser, ("builtin", "cyvcf2"), "parser")
    flag = _flag(filter, "filter")
    vcf = input_path(vcf_path)
    popinfo = input_path(popinfo_path)
    try:
        if mode == "cyvcf2":
            require("cyvcf2", _DADI_EXTRA)
            data = dadi.Misc.make_data_dict_vcf_cyvcf2(str(vcf), str(popinfo), filter=flag)
        else:
            data = dadi.Misc.make_data_dict_vcf(str(vcf), str(popinfo), filter=flag)
    except (ValueError, KeyError, IndexError, OSError) as exc:
        raise ConfigurationError(f"dadi could not read {vcf.name} as VCF: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_vcf_data",
        "parser": mode,
        "source": str(vcf),
        "popinfo": str(popinfo),
        "populations": _populations(data),
        "num_snps": len(data),
        "provenance": _provenance("dadi", f"dadi VCF reader ({mode}); biallelic SNPs"),
    }


def bootstrap_snp_spectra(
    path: str | Path,
    *,
    pop_ids: Sequence[str],
    projections: Sequence[int],
    chunk_size: int,
    num_boot: int,
    seed: int = 0,
    polarized: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bootstrap frequency spectra from chunks of an SNP file, using dadi's bootstrap helpers.

    The SNP IDs must end in the base-pair position (for example ``chr1_1500``). The file is cut
    into chunks of ``chunk_size`` base pairs, and each bootstrap draws chunks with replacement.
    ``projections`` gives the projected sample size of each population. The result is the mean
    and standard deviation of the bootstrap spectra. ``seed`` makes the draws reproducible.
    """

    dadi, np = _dadi_and_numpy()
    source = input_path(path)
    data = _load_snp(dadi, source)
    selected = _pop_ids(pop_ids, _populations(data))
    if len(selected) > 3:
        raise ConfigurationError("bootstrap_snp_spectra supports one to three populations.")
    proj = [_integer(item, "projection", minimum=2, maximum=_MAX_2D_SAMPLE) for item in projections]
    if len(proj) != len(selected):
        raise ConfigurationError("projections must give one sample size per population.")
    size = _integer(chunk_size, "chunk_size", minimum=1)
    boots = _integer(num_boot, "num_boot", minimum=1, maximum=_MAX_BOOTSTRAPS)
    seed_value = _integer(seed, "seed", minimum=0)
    polar = _flag(polarized, "polarized")
    try:
        fragments = dadi.Misc.fragment_data_dict(data, size)
    except (ValueError, KeyError, IndexError) as exc:
        raise ConfigurationError(
            f"SNP IDs must end in _<position> to form chunks of {size} bp: {exc}"
        ) from exc
    state = random.getstate()
    try:
        random.seed(seed_value)
        spectra = dadi.Misc.bootstraps_from_dd_chunks(
            fragments, boots, selected, proj, mask_corners=True, polarized=polar
        )
    except (ValueError, KeyError, IndexError) as exc:
        raise ConfigurationError(f"dadi could not bootstrap {source.name}: {exc}") from exc
    finally:
        random.setstate(state)
    stack = np.array(
        [np.ma.filled(np.ma.asarray(spectrum, dtype=float), np.nan) for spectrum in spectra]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(stack, axis=0)
        spread = np.nanstd(stack, axis=0)
    mean_masked = np.ma.masked_invalid(mean)
    spread_masked = np.ma.masked_invalid(spread)
    return {
        "status": "completed",
        "operation": "bootstrap_snp_spectra",
        "source": str(source),
        "pop_ids": selected,
        "projections": proj,
        "chunk_size": size,
        "num_chunks": len(fragments),
        "num_boot": boots,
        "seed": seed_value,
        "mean_spectrum": _nested(np, mean_masked),
        "std_spectrum": _nested(np, spread_masked),
        "provenance": _provenance("dadi", "bootstrap of chunked SNP data; masks as dadi"),
    }


def pool_sfs_populations(
    spectrum: Sequence[Any],
    *,
    populations: Sequence[int] = (0, 1),
    verbose: bool = True,
) -> dict[str, Any]:
    """Pool two populations of a two- or three-population spectrum, using dadi's ``combine_pops()``.

    The two populations named by ``populations`` (0-based axes) are summed into one axis, which
    becomes the first axis of the result. For three populations the pair must be (0, 1), (0, 2)
    or (1, 2), as dadi supports.
    """

    dadi, np = _dadi_and_numpy()
    counts = _counts_array(np, spectrum, "spectrum", dims=(2, 3), nonnegative=False)
    pair = _pair(populations, counts.ndim)
    if counts.ndim == 3 and pair not in ((0, 1), (0, 2), (1, 2)):
        raise ConfigurationError(
            "for three populations, populations must be (0, 1), (0, 2) or (1, 2)."
        )
    fs = dadi.Spectrum(counts)
    ns = _spectrum_sizes(fs)
    pooled = dadi.Misc.combine_pops(fs, idx=list(pair))
    result_sizes = [ns[pair[0]] + ns[pair[1]]] + [
        ns[axis] for axis in range(len(ns)) if axis not in pair
    ]
    return {
        "status": "completed",
        "operation": "pool_sfs_populations",
        "populations": list(pair),
        "sample_sizes": result_sizes,
        "spectrum": _nested(np, dadi.Spectrum(np.asarray(pooled, dtype=float))),
        "provenance": _provenance("dadi", "dadi combine_pops: pooled populations are summed"),
    }


def _pair(value: object, ndim: int) -> tuple[int, int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError("populations must be two distinct 0-based population indices.")
    first = _integer(value[0], "population index", minimum=0, maximum=ndim - 1)
    second = _integer(value[1], "population index", minimum=0, maximum=ndim - 1)
    if first == second:
        raise ConfigurationError("populations must be two distinct indices.")
    return (first, second)


def read_sfs_file(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read a dadi array or spectrum text file, using dadi's ``Numerics.array_from_file()``.

    The first non-comment line gives the shape and the following lines give the values. Lines
    starting with ``#`` are returned as comments.
    """

    dadi, np = _dadi_and_numpy()
    source = input_path(path)
    try:
        with source.open(encoding="utf8") as handle:
            array, comments = dadi.Numerics.array_from_file(handle, return_comments=True)
    except (ValueError, OSError, IndexError) as exc:
        raise ConfigurationError(f"dadi could not read {source.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_sfs_file",
        "source": str(source),
        "shape": [int(n) for n in np.asarray(array).shape],
        "spectrum": _nested(np, np.ma.masked_array(np.asarray(array, dtype=float))),
        "comments": [str(line) for line in comments],
        "provenance": _provenance("dadi", "dadi array_from_file text format"),
    }


def write_sfs_file(
    spectrum: Sequence[Any],
    destination: str | Path,
    *,
    precision: int = 16,
    comment_lines: Sequence[str] = (),
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a spectrum or array to a dadi text file, using dadi's ``Numerics.array_to_file()``.

    The file starts with the shape, then the values. ``comment_lines`` are written after ``# ``
    at the top of the file. An existing file is never overwritten.
    """

    dadi, np = _dadi_and_numpy()
    array = _float_array(np, spectrum, "spectrum")
    digits = _integer(precision, "precision", minimum=1, maximum=17)
    comments = [_text(line, "comment line") for line in comment_lines]
    target = _new_file(destination, "destination")
    with target.open("w", encoding="utf8") as handle:
        dadi.Numerics.array_to_file(array, handle, precision=digits, comment_lines=comments)
    return {
        "status": "completed",
        "operation": "write_sfs_file",
        "destination": str(target),
        "bytes": target.stat().st_size,
        "shape": [int(n) for n in array.shape],
        "num_values": int(array.size),
        "provenance": _provenance("dadi", "dadi array_to_file text format"),
    }


def _float_array(np: Any, value: object, name: str) -> Any:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular list of numbers.") from exc
    if array.ndim < 1 or array.size == 0 or not np.isfinite(array).all():
        raise ConfigurationError(f"{name} must be a non-empty array of finite numbers.")
    return array


def apply_ancestral_misidentification(
    spectrum: Sequence[Any],
    p_misid: float,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply ancestral-state misidentification to a spectrum, using dadi's anc_state_misid helper.

    The result is ``(1 - p) * spectrum + p * reversed(spectrum)``, where ``p`` is the fraction of
    sites whose ancestral state is taken to be wrong.
    """

    dadi, np = _dadi_and_numpy()
    counts = _counts_array(np, spectrum, "spectrum", nonnegative=False)
    prob = _number(p_misid, "p_misid", minimum=0.0)
    if prob > 1.0:
        raise ConfigurationError("p_misid must be between 0 and 1.")
    fs = dadi.Spectrum(counts)
    result = dadi.Numerics.apply_anc_state_misid(fs, prob)
    return {
        "status": "completed",
        "operation": "apply_ancestral_misidentification",
        "p_misid": prob,
        "sample_sizes": _spectrum_sizes(fs),
        "spectrum": _nested(np, result),
        "provenance": _provenance("dadi", "(1 - p) * fs + p * reversed fs"),
    }


def _demographic_expected(
    np: Any,
    family_module: Any,
    options: Sequence[str],
    model: object,
    params: object,
    ns: tuple[int, ...],
    points: int,
    operation: str,
) -> dict[str, Any]:
    name = _choice(model, options, "model")
    function = getattr(family_module, name, None)
    if function is None:
        raise ConfigurationError(f"Installed dadi does not provide model {name!r}.")
    values = _numbers(params, "params")
    spectrum = _run_model(function, values, ns, points)
    return {
        "status": "completed",
        "operation": operation,
        "model": name,
        "params": values,
        "sample_sizes": list(ns),
        "grid_points": points,
        "spectrum": _nested(np, spectrum),
        "provenance": _provenance("dadi", _MODEL_CONVENTION + "; expected under the model"),
    }


def expected_sfs_1d(
    model: str,
    params: Sequence[float],
    sample_size: int,
    *,
    grid_points: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected site frequency spectrum for a one-population demographic model, using dadi.

    ``model`` picks a model from the one-population family: ``two_epoch`` (the default
    representative), ``growth``, ``bottlegrowth``, ``three_epoch`` or ``three_epoch_inbreeding``.
    ``params`` are that model's parameters, in dadi's order. ``grid_points`` is the integration
    grid; the result approaches the exact spectrum as it grows.
    """

    dadi, np = _dadi_and_numpy()
    n = _integer(sample_size, "sample_size", minimum=2, maximum=_MAX_1D_SAMPLE)
    points = _integer(grid_points, "grid_points", minimum=50, maximum=10_000)
    return _demographic_expected(
        np,
        dadi.Demographics1D,
        _DEMOGRAPHICS_1D,
        model,
        params,
        (n,),
        points,
        "expected_sfs_1d",
    )


def expected_sfs_2d(
    model: str,
    params: Sequence[float],
    sample_sizes: Sequence[int],
    *,
    grid_points: int = 60,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected joint spectrum for a two-population demographic model, using dadi.

    ``model`` picks one of dadi's two-population models, for example ``split_mig`` (the default
    representative), ``asym_mig``, ``anc_sym_mig`` or ``IM``. The other names of the family are
    accepted too, and each has the same ``(params, ns, pts)`` interface.
    """

    dadi, np = _dadi_and_numpy()
    ns = _sizes(sample_sizes, 2, _MAX_2D_SAMPLE)
    points = _integer(grid_points, "grid_points", minimum=20, maximum=400)
    return _demographic_expected(
        np,
        dadi.Demographics2D,
        _DEMOGRAPHICS_2D,
        model,
        params,
        ns,
        points,
        "expected_sfs_2d",
    )


def expected_sfs_3d(
    model: str,
    params: Sequence[float],
    sample_sizes: Sequence[int],
    *,
    grid_points: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected joint spectrum for a three-population demographic model, using dadi.

    ``model`` picks one of dadi's three-population models; the default representative is
    ``split_nomig``, and ``out_of_africa`` is the Gutenkunst et al. (2009) model. The grid is
    kept small because the three-population integration grows as the cube of ``grid_points``.
    """

    dadi, np = _dadi_and_numpy()
    ns = _sizes(sample_sizes, 3, _MAX_3D_SAMPLE)
    points = _integer(grid_points, "grid_points", minimum=20, maximum=60)
    return _demographic_expected(
        np,
        dadi.Demographics3D,
        _DEMOGRAPHICS_3D,
        model,
        params,
        ns,
        points,
        "expected_sfs_3d",
    )


def expected_sfs_selection(
    model: str,
    params: Sequence[float],
    sample_sizes: Sequence[int],
    *,
    grid_points: int = 40,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected spectrum for a demographic model with purifying or positive selection, using dadi.

    These are dadi's DFE models, which add selection coefficients (``gamma``) to a demographic
    history. ``model`` picks the model, for example ``split_mig_sel`` (the default
    representative); the ``*_single_gamma`` variants share one selection coefficient across
    populations. ``sample_sizes`` has one or two entries, as the chosen model requires.
    """

    np = require("numpy", _DADI_EXTRA)
    if isinstance(sample_sizes, (str, bytes)) or not isinstance(sample_sizes, Sequence):
        raise ConfigurationError("sample_sizes must list one or two integer sample sizes.")
    if not 1 <= len(sample_sizes) <= 2:
        raise ConfigurationError("sample_sizes must list one or two integer sample sizes.")
    ns = tuple(
        _integer(item, "sample size", minimum=2, maximum=_MAX_2D_SAMPLE) for item in sample_sizes
    )
    points = _integer(grid_points, "grid_points", minimum=20, maximum=400)
    name = _choice(model, _DEMOGRAPHICS_SELECTION, "model")
    selection = importlib.import_module("dadi.DFE.DemogSelModels")
    function = getattr(selection, name, None)
    if function is None:
        raise ConfigurationError(f"Installed dadi does not provide selection model {name!r}.")
    values = _numbers(params, "params")
    spectrum = _run_model(function, values, ns, points)
    return {
        "status": "completed",
        "operation": "expected_sfs_selection",
        "model": name,
        "params": values,
        "sample_sizes": list(ns),
        "grid_points": points,
        "spectrum": _nested(np, spectrum),
        "provenance": _provenance("dadi", _MODEL_CONVENTION + "; demography with selection"),
    }


def expected_sfs_from_demes(
    graph: str | Path,
    sampled_demes: Sequence[str],
    sample_sizes: Sequence[int],
    *,
    grid_points: int = 40,
    theta: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected spectrum for a demes demographic graph file, using dadi's ``Demes.SFS()``.

    ``graph`` is a demes YAML file (read with the ``demes`` package, which must be installed).
    ``sampled_demes`` names one to three demes of the graph, and ``sample_sizes`` gives the
    number of chromosomes sampled from each. ``theta`` scales the spectrum.
    """

    dadi, np = _dadi_and_numpy()
    path = input_path(graph)
    if isinstance(sampled_demes, (str, bytes)) or not isinstance(sampled_demes, Sequence):
        raise ConfigurationError("sampled_demes must list one to three deme names.")
    if not 1 <= len(sampled_demes) <= 3:
        raise ConfigurationError("sampled_demes must list one to three deme names.")
    names = [_text(name, "deme name") for name in sampled_demes]
    if len(set(names)) != len(names):
        raise ConfigurationError("sampled_demes must not repeat a deme.")
    maximum = {1: _MAX_1D_SAMPLE, 2: _MAX_2D_SAMPLE, 3: _MAX_3D_SAMPLE}[len(names)]
    ns = _sizes(sample_sizes, len(names), maximum)
    points = _integer(
        grid_points, "grid_points", minimum=20, maximum=60 if len(names) == 3 else 400
    )
    scale = _number(theta, "theta", positive=True)
    demes = _demes_module()
    try:
        demography = demes.load(str(path))
    except (ValueError, KeyError, TypeError, OSError) as exc:
        raise ConfigurationError(f"demes could not load {path.name}: {exc}") from exc
    try:
        spectrum = dadi.Demes.SFS(demography, names, ns, points, theta=scale)
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        raise ConfigurationError(f"dadi could not compute the demes spectrum: {exc}") from exc
    return {
        "status": "completed",
        "operation": "expected_sfs_from_demes",
        "source": str(path),
        "demes": names,
        "sample_sizes": list(ns),
        "grid_points": points,
        "theta": scale,
        "spectrum": _nested(np, spectrum),
        "provenance": _provenance("dadi", _MODEL_CONVENTION + "; demes graph via dadi Demes.SFS"),
    }


def _demes_module() -> Any:
    try:
        return importlib.import_module("demes")
    except ImportError as exc:
        raise BackendUnavailableError(
            f"demes could not be imported: {exc}",
            hint='Install the demes package with: python -m pip install "demes"',
        ) from exc


def fit_demographic_model(
    data: Sequence[Any],
    *,
    model: str,
    p0: Sequence[float],
    lower_bound: Sequence[float] | None = None,
    upper_bound: Sequence[float] | None = None,
    method: str = "bfgs",
    multinom: bool = True,
    grid_points: int = 40,
    maxiter: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a demographic model to an observed spectrum by maximum likelihood, using dadi.

    ``data`` is the observed spectrum (one or two populations). ``model`` is a model from the
    one-population or two-population family, chosen by the dimension of ``data``. ``p0`` are the
    starting parameters and must be positive. ``method`` picks the optimizer: ``bfgs`` (BFGS on
    log parameters), ``lbfgsb``, ``nelder-mead`` or ``powell``. With ``multinom=True`` dadi
    rescales theta, so ``p0`` does not include it; the log likelihood is reported at the fitted
    parameters.
    """

    dadi, np = _dadi_and_numpy()
    counts = _counts_array(np, data, "data", dims=(1, 2))
    observed = dadi.Spectrum(counts)
    ns = tuple(_spectrum_sizes(observed))
    family, options = (
        (dadi.Demographics1D, _DEMOGRAPHICS_1D)
        if counts.ndim == 1
        else (dadi.Demographics2D, _DEMOGRAPHICS_2D)
    )
    function = _model_function(family, model, options)
    optimizer_name = _OPTIMIZERS[_choice(method, _FIT_OPTIMIZER_METHODS, "method")]
    optimizer = getattr(dadi.Inference, optimizer_name)
    start = _numbers(p0, "p0")
    if any(value <= 0 for value in start):
        raise ConfigurationError("p0 must contain only positive values.")
    lower = _bounds(lower_bound, len(start), "lower_bound")
    upper = _bounds(upper_bound, len(start), "upper_bound")
    if (
        lower is not None
        and upper is not None
        and any(lo >= hi for lo, hi in zip(lower, upper, strict=True))
    ):
        raise ConfigurationError("each lower_bound must be below its upper_bound.")
    flag = _flag(multinom, "multinom")
    points = _integer(grid_points, "grid_points", minimum=20, maximum=400)
    iterations = _integer(maxiter, "maxiter", minimum=1, maximum=_MAX_MAXITER)
    try:
        fitted = optimizer(
            start,
            observed,
            function,
            points,
            lower_bound=lower,
            upper_bound=upper,
            verbose=0,
            multinom=flag,
            maxiter=iterations,
        )
        model_spectrum = function(tuple(float(value) for value in fitted), ns, points)
        if flag:
            log_likelihood = dadi.Inference.ll_multinom(model_spectrum, observed)
        else:
            log_likelihood = dadi.Inference.ll(model_spectrum, observed)
    except (ValueError, TypeError, IndexError, KeyError) as exc:
        if "iprint" in str(exc):
            raise ConfigurationError(
                "dadi's L-BFGS-B optimizer passes 'iprint' to SciPy, which newer SciPy releases "
                "(1.18 and later) no longer accept; use method='bfgs', 'nelder-mead' or 'powell'."
            ) from exc
        raise ConfigurationError(f"dadi could not fit {model!r} to the data: {exc}") from exc
    return {
        "status": "completed",
        "operation": "fit_demographic_model",
        "model": model,
        "method": method,
        "multinom": flag,
        "sample_sizes": list(ns),
        "grid_points": points,
        "maxiter": iterations,
        "params": [float(value) for value in np.atleast_1d(fitted)],
        "log_likelihood": _finite(log_likelihood),
        "provenance": _provenance("dadi", f"maximum likelihood; dadi {optimizer_name}"),
    }


def sfs_residuals(
    model: Sequence[Any],
    data: Sequence[Any],
    *,
    residual: str = "anscombe",
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-entry residuals between a model spectrum and observed counts, using dadi.

    ``residual="anscombe"`` gives dadi's Anscombe Poisson residuals; ``"poisson"`` gives the
    linear Poisson residuals ``(model - data) / sqrt(model)``. Both arrays must have the same
    shape. Masked corners are None, and ``sum_of_squares`` covers the unmasked entries.
    """

    dadi, np = _dadi_and_numpy()
    kind = _choice(residual, tuple(_RESIDUALS), "residual")
    model_array = _counts_array(np, model, "model", nonnegative=False)
    data_array = _counts_array(np, data, "data")
    if model_array.shape != data_array.shape:
        raise ConfigurationError("model and data must have the same shape.")
    residuals = getattr(dadi.Inference, _RESIDUALS[kind])(
        dadi.Spectrum(model_array), dadi.Spectrum(data_array)
    )
    masked = np.ma.masked_invalid(np.ma.asarray(residuals, dtype=float))
    values = np.ma.compressed(masked)
    return {
        "status": "completed",
        "operation": "sfs_residuals",
        "residual": kind,
        "shape": [int(n) for n in model_array.shape],
        "residuals": _nested(np, masked),
        "num_entries": int(values.size),
        "sum_of_squares": _finite(float(np.sum(values**2))) if values.size else 0.0,
        "provenance": _provenance("dadi", f"dadi {_RESIDUALS[kind]}"),
    }


def parameter_uncertainty_fim(
    data: Sequence[Any],
    *,
    model: str,
    p0: Sequence[float],
    grid_points: Sequence[int] = (20, 30, 40),
    log: bool = False,
    multinom: bool = True,
    eps: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parameter standard errors from the Fisher information matrix, using dadi's Godambe module.

    The model is extrapolated over ``grid_points`` (three grid sizes by default), as dadi
    recommends, and its Hessian of the log likelihood gives the Fisher information. With
    ``multinom=True`` theta is estimated too, and its uncertainty is the last entry. ``log=True``
    gives uncertainties on the log scale. ``eps`` is the relative step of the finite differences.
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
    start = _numbers(p0, "p0")
    if any(value <= 0 for value in start):
        raise ConfigurationError("p0 must contain only positive values.")
    grids = [_integer(item, "grid_points entry", minimum=20, maximum=400) for item in grid_points]
    if not grids or len(grids) > 5:
        raise ConfigurationError("grid_points must list one to five grid sizes.")
    flag_log = _flag(log, "log")
    flag_multi = _flag(multinom, "multinom")
    step = _number(eps, "eps", positive=True)
    if step >= 1.0:
        raise ConfigurationError("eps must be below 1.")
    try:
        extrapolated = dadi.Numerics.make_extrap_func(function)
        uncertainties = dadi.Godambe.FIM_uncert(
            extrapolated,
            grids,
            start,
            observed,
            log=flag_log,
            multinom=flag_multi,
            eps=step,
        )
    except (ValueError, TypeError, IndexError, KeyError) as exc:
        raise ConfigurationError(f"dadi could not compute the Fisher information: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parameter_uncertainty_fim",
        "model": model,
        "params": start,
        "grid_points": grids,
        "log": flag_log,
        "multinom": flag_multi,
        "uncertainties": [_finite(float(value)) for value in np.atleast_1d(uncertainties)],
        "provenance": _provenance("dadi", "Fisher information from a numerical Hessian (Godambe)"),
    }


def annotate_snps_from_annovar(
    path: str | Path,
    annovar_path: str | Path,
    *,
    variant_type: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the SNPs of one annotated variant type, using dadi and an ANNOVAR output file.

    ``annovar_path`` is an ``.exonic_variant_function`` file from ANNOVAR. A SNP is kept when its
    variant type text (for example ``nonsynonymous SNV``) contains ``variant_type``. The SNP IDs
    of the kept sites are returned in sorted order.
    """

    dadi = require("dadi", _DADI_EXTRA)
    source = input_path(path)
    annotation = input_path(annovar_path)
    label = _text(variant_type, "variant_type")
    data = _load_snp(dadi, source)
    try:
        kept = dadi.Misc.annotate_from_annovar(data, str(annotation), label)
    except (ValueError, IndexError, OSError) as exc:
        raise ConfigurationError(f"dadi could not read {annotation.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "annotate_snps_from_annovar",
        "source": str(source),
        "annovar": str(annotation),
        "variant_type": label,
        "num_snps_total": len(data),
        "num_snps_kept": len(kept),
        "kept_snp_ids": sorted(kept),
        "provenance": _provenance("dadi", "ANNOVAR .exonic_variant_function variant type filter"),
    }


def trinucleotide_frequency_array(
    frequencies: Mapping[str, float],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert trinucleotide frequencies to dadi's 64-entry array, using dadi's helper.

    ``frequencies`` maps each of the 64 trinucleotides (over A, C, G and T) to a non-negative
    frequency. dadi orders them with the base sequence ``C, G, T, A`` for each of the three
    positions, the first position varying slowest; ``trinucleotides`` gives that order.
    """

    dadi, np = _dadi_and_numpy()
    if not isinstance(frequencies, Mapping):
        raise ConfigurationError("frequencies must map each trinucleotide to a frequency.")
    order = [
        left + center + right
        for left in _TRI_ORDER_BASES
        for center in _TRI_ORDER_BASES
        for right in _TRI_ORDER_BASES
    ]
    keyset = set("ACGT")
    values: dict[str, float] = {}
    for key in sorted(frequencies):
        if not isinstance(key, str) or len(key) != 3 or set(key) - keyset:
            raise ConfigurationError(f"{key!r} is not a trinucleotide over A, C, G and T.")
        values[key] = _number(frequencies[key], f"frequency of {key}", minimum=0.0)
    missing = sorted(set(order) - set(values))
    if missing:
        raise ConfigurationError(f"frequencies is missing trinucleotides {missing[:5]}.")
    array = dadi.Misc.tri_freq_dict_to_array({key: values[key] for key in order})
    return {
        "status": "completed",
        "operation": "trinucleotide_frequency_array",
        "trinucleotides": order,
        "values": [float(value) for value in np.asarray(array, dtype=float)],
        "provenance": _provenance("dadi", "dadi tri_freq_dict_to_array order C, G, T, A"),
    }


def _dadi_and_numpy() -> tuple[Any, Any]:
    dadi = require("dadi", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    return dadi, np


__all__ = [
    "annotate_snps_from_annovar",
    "apply_ancestral_misidentification",
    "bootstrap_snp_spectra",
    "expected_sfs_1d",
    "expected_sfs_2d",
    "expected_sfs_3d",
    "expected_sfs_from_demes",
    "expected_sfs_selection",
    "fit_demographic_model",
    "parameter_uncertainty_fim",
    "pool_sfs_populations",
    "read_sfs_file",
    "read_snp_data",
    "read_vcf_data",
    "sfs_residuals",
    "snp_configuration_counts",
    "trinucleotide_frequency_array",
    "write_sfs_file",
]

_register_public_verbose(globals(), tuple(__all__))
