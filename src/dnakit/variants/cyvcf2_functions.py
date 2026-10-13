"""Optional cyvcf2 adapters for genotype heterozygosity, relatedness and correlation.

- ``summarize_sample_heterozygosity()`` wraps ``cyvcf2.par_het``: per-sample heterozygous fraction
  and allele-balance percentiles at explicit sites.
- ``summarize_pairwise_relatedness()`` wraps ``cyvcf2.par_relatedness``: pairwise IBS counts and a
  relatedness coefficient at explicit sites.
- ``pearson_genotype_correlation()`` wraps ``cyvcf2.r_unphased``: Pearson correlation between two
  dosage vectors, with missing calls skipped pairwise.

The VCF must be bgzip-compressed with a tabix or CSI index, because the upstream functions query
each site by region. Sites are always explicit, so no upstream thinning rule chooses them. Upstream
still skips some sites (call rate below 0.8, or allele frequency below 0.01 in the heterozygosity
path); the ``sites_used`` field lists the sites that were kept.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``bio-cyvcf2`` install hint.
"""

from __future__ import annotations

import importlib.metadata
import math
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-cyvcf2"
_INDEX_SUFFIXES = (".tbi", ".csi")
_ALLELES = frozenset("ACGTN")
_MAX_SITES = 100_000
_MAX_SAMPLES = 10_000
_MAX_NCPUS = 8
_MAX_VALUES = 1_000_000
_MAX_POSITION = 2_147_483_647
_UNKNOWN_CODE = 3  # cyvcf2's internal code for a missing call


def _cyvcf2() -> Any:
    return require("cyvcf2", _EXTRA)


def _provenance(convention: str) -> dict[str, str]:
    return {
        "package": "cyvcf2",
        "version": importlib.metadata.version("cyvcf2"),
        "convention": convention,
    }


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _sequence(value: object, name: str, *, limit: int) -> list[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    if len(value) > limit:
        raise ConfigurationError(f"{name} exceeds the limit of {limit} items.")
    return list(value)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n"):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value.strip()


def _allele(value: object, name: str) -> str:
    text = _text(value, name).upper()
    if not set(text) <= _ALLELES:
        raise ConfigurationError(f"{name} must contain only A, C, G, T or N.")
    return text


def _indexed_vcf(value: object) -> Path:
    if isinstance(value, Path):
        value = str(value)
    path = input_path(_text(value, "source"))
    if path.suffix != ".gz":
        raise ConfigurationError(f"{path.name} must be bgzip-compressed (.gz) and indexed.")
    if not any(Path(str(path) + suffix).is_file() for suffix in _INDEX_SUFFIXES):
        raise ConfigurationError(f"{path.name} needs a tabix (.tbi) or CSI (.csi) index.")
    return path


def _sample_names(path: Path, value: object) -> list[str]:
    names = [_text(item, "samples") for item in _sequence(value, "samples", limit=_MAX_SAMPLES)]
    if any("," in name for name in names) or len(set(names)) != len(names):
        raise ConfigurationError("samples must be unique.")
    header = _cyvcf2().VCF(str(path)).samples
    known = set(header)
    missing = [name for name in names if name not in known]
    if missing:
        raise ConfigurationError(f"unknown sample names: {missing}")
    return names


def _sites(value: object) -> list[list[Any]]:
    result: list[list[Any]] = []
    for item in _sequence(value, "sites", limit=_MAX_SITES):
        if not isinstance(item, (list, tuple)) or len(item) not in (2, 4):
            raise ConfigurationError("each site must be (chrom, pos) or (chrom, pos, ref, alt).")
        chrom = _text(item[0], "site chrom")
        if "," in chrom:
            raise ConfigurationError("site chromosome names must not contain commas.")
        if chrom in ("X", "chrX"):
            raise ConfigurationError("upstream heterozygosity stops at chromosome X; remove it.")
        pos = _integer(item[1], "site pos", minimum=1, maximum=_MAX_POSITION)
        if len(item) == 2:
            result.append([chrom, pos])
        else:
            result.append([chrom, pos, _allele(item[2], "site ref"), _allele(item[3], "site alt")])
    return result


def _plain(value: object) -> Any:
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return float(value) if math.isfinite(value) else None
    return str(value)


def _messages(caught: list[warnings.WarningMessage]) -> list[str]:
    # Upstream never closes its multiprocessing pool; that resource notice is not analysis output.
    texts = {str(item.message) for item in caught}
    return sorted(text for text in texts if "multiprocessing pool" not in text)


def summarize_sample_heterozygosity(
    source: str | Path,
    samples: Sequence[str],
    sites: Sequence[Sequence[object]],
    *,
    ncpus: int = 1,
    min_depth: int = 8,
    percentiles: Sequence[int] = (10, 90),
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-sample heterozygosity and allele-balance percentiles at explicit sites.

    ``sites`` holds ``(chrom, pos)`` or ``(chrom, pos, ref, alt)`` tuples. Depths come from the
    FORMAT AD field, and only depths above ``min_depth`` enter the percentiles.
    """

    path = _indexed_vcf(source)
    names = _sample_names(path, samples)
    site_list = _sites(sites)
    workers = _integer(ncpus, "ncpus", minimum=1, maximum=_MAX_NCPUS)
    depth = _integer(min_depth, "min_depth", minimum=0, maximum=100_000)
    if isinstance(percentiles, (str, bytes)) or len(percentiles) != 2:
        raise ConfigurationError("percentiles must be a pair of integers such as (10, 90).")
    low, high = (_integer(item, "percentile", minimum=0, maximum=100) for item in percentiles)
    if low >= high:
        raise ConfigurationError("percentiles must be two increasing integers in [0, 100].")
    upstream = _cyvcf2()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ranges, used, _genotypes = upstream.par_het(
            str(path), names, workers, site_list, min_depth=depth, percentiles=(low, high)
        )
    samples_out = {
        name: {key: _plain(value) for key, value in sorted(ranges[name].items())}
        for name in sorted(ranges)
    }
    return {
        "status": "completed",
        "operation": "summarize_sample_heterozygosity",
        "samples": samples_out,
        "sites_requested": len(site_list),
        "sites_used": sorted(str(item) for item in used),
        "warnings": _messages(caught),
        "provenance": _provenance(
            "het_ratio = het_count / sites_used; percentiles of alt-allele fraction "
            "over heterozygous calls with depth above min_depth"
        ),
    }


def summarize_pairwise_relatedness(
    source: str | Path,
    samples: Sequence[str],
    sites: Sequence[Sequence[object]],
    *,
    ncpus: int = 1,
    min_depth: int = 5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise IBS counts and relatedness at explicit sites, as in cyvcf2's relatedness code.

    Each unordered sample pair appears once, sorted by sample name. ``relatedness`` is the
    upstream coefficient (about 0.5 for parent-offspring, 1 for identical samples in typical
    cohorts). A sample with no heterozygous sites produces a warning and an unreliable value.
    """

    path = _indexed_vcf(source)
    names = _sample_names(path, samples)
    site_list = _sites(sites)
    workers = _integer(ncpus, "ncpus", minimum=1, maximum=_MAX_NCPUS)
    depth = _integer(min_depth, "min_depth", minimum=0, maximum=100_000)
    upstream = _cyvcf2()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        arrays = upstream.par_relatedness(str(path), names, workers, site_list, min_depth=depth)
    pairs = [
        {
            "sample_a": str(arrays["sample_a"][index]),
            "sample_b": str(arrays["sample_b"][index]),
            "relatedness": _plain(arrays["rel"][index]),
            "ibs0": _plain(arrays["ibs0"][index]),
            "ibs2": _plain(arrays["ibs2"][index]),
            "informative_sites": _plain(arrays["n"][index]),
            "hets_a": _plain(arrays["hets_a"][index]),
            "hets_b": _plain(arrays["hets_b"][index]),
            "shared_hets": _plain(arrays["shared_hets"][index]),
        }
        for index in range(len(arrays["sample_a"]))
    ]
    pairs.sort(key=lambda item: (item["sample_a"], item["sample_b"]))
    return {
        "status": "completed",
        "operation": "summarize_pairwise_relatedness",
        "pairs": pairs,
        "sites_requested": len(site_list),
        "warnings": _messages(caught),
        "provenance": _provenance(
            "relatedness = (ibs0 - 2 * ibs2) / min(hets_a, hets_b), as computed by cyvcf2"
        ),
    }


def pearson_genotype_correlation(
    first: Sequence[int | None],
    second: Sequence[int | None],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pearson correlation between two genotype dosage vectors (0, 1, 2; None is missing).

    Pairs where either call is missing are skipped. The correlation is None when fewer than two
    usable pairs remain or when either vector is constant over those pairs.
    """

    left = _sequence(list(first), "first", limit=_MAX_VALUES)
    right = _sequence(list(second), "second", limit=_MAX_VALUES)
    if len(left) != len(right):
        raise ConfigurationError("first and second must have the same length.")
    for value in (*left, *right):
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value not in (0, 1, 2):
            raise ConfigurationError("genotype dosages must be 0, 1, 2 or None.")
    usable = [(a, b) for a, b in zip(left, right, strict=True) if a is not None and b is not None]
    constant = len({a for a, _ in usable}) < 2 or len({b for _, b in usable}) < 2
    value: float | None = None
    if len(usable) >= 2 and not constant:
        np = require("numpy", "bio-analysis")
        upstream = _cyvcf2()
        a_codes = np.array(
            [_UNKNOWN_CODE if item is None else item for item in left], dtype=np.int32
        )
        b_codes = np.array(
            [_UNKNOWN_CODE if item is None else item for item in right], dtype=np.int32
        )
        # The upstream frequency argument is unused; 0.0 is passed to keep the signature.
        raw = upstream.r_unphased(a_codes, b_codes, 0.0, len(left))
        value = _plain(raw)
    return {
        "status": "completed",
        "operation": "pearson_genotype_correlation",
        "correlation": value,
        "usable_pairs": len(usable),
        "missing_pairs": len(left) - len(usable),
        "provenance": _provenance(
            "Pearson r over pairs where both calls are present; the upstream frequency argument "
            "is ignored"
        ),
    }


__all__ = [
    "pearson_genotype_correlation",
    "summarize_pairwise_relatedness",
    "summarize_sample_heterozygosity",
]

_register_public_verbose(globals(), tuple(__all__))
