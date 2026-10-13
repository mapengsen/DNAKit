"""Optional-backend VCF readers and scikit-allel population statistics.

The functions in this module are explicit adapters over optional packages:

- ``read_vcf()`` and ``fetch_vcf_region()`` read VCF files with ``cyvcf2``,
  ``pysam`` or ``vcfpy``. Every backend returns the same record shape, so the
  output can be passed directly to the dependency-free ``calculate_*`` helpers
  in :mod:`dnakit.variants.statistics`.
- ``calculate_*`` population statistics use ``scikit-allel`` and read the same
  record shape, with genotype calls such as ``"0/1"`` or ``"1|1"``.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.variants.statistics import _base, _records

_IO_BACKENDS = ("cyvcf2", "pysam", "vcfpy")
_REGION_BACKENDS = ("cyvcf2", "pysam")
_EXTRAS = {"cyvcf2": "bio-cyvcf2", "pysam": "bio-reads", "vcfpy": "bio-vcfpy"}
_NUMPY_EXTRA = "bio-analysis"
_SCIKIT_ALLEL_EXTRA = "bio-scikit-allel"
_DEFAULT_MAX_RECORDS = 100_000
_MAX_RECORDS_LIMIT = 1_000_000
_REGION_PATTERN = re.compile(r"^(?P<chrom>[^:\s]+)(?::(?P<start>[\d,]+)(?:-(?P<end>[\d,]+))?)?$")
_IGNORED_FILTERS = frozenset({"", ".", "PASS"})


@dataclass(frozen=True)
class _Panel:
    """Genotype calls for the sites that survived the requested filters."""

    samples: list[str]
    chromosome: str | None
    positions: Any
    genotypes: Any
    skipped: dict[str, int]


def _choose(value: object, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"backend must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _limit(value: object) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 1 <= value <= _MAX_RECORDS_LIMIT
    ):
        raise ConfigurationError(f"max_records must be an integer in [1, {_MAX_RECORDS_LIMIT}].")
    return value


def _position(text: str) -> int:
    value = int(text.replace(",", ""))
    if value < 1:
        raise ConfigurationError("region positions are 1-based and must be at least 1.")
    return value


def _parse_region(region: object) -> tuple[str, int | None, int | None]:
    match = _REGION_PATTERN.match(region.strip()) if isinstance(region, str) else None
    if match is None:
        raise ConfigurationError(
            "region must look like 'chr1', 'chr1:1000' or 'chr1:1000-2000' (1-based, inclusive)."
        )
    start = _position(match["start"]) if match["start"] else None
    end = _position(match["end"]) if match["end"] else start
    if start is not None and end is not None and end < start:
        raise ConfigurationError("region end must be greater than or equal to its start.")
    return match["chrom"], start, end


def _json_value(value: object) -> object:
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        # VCF numbers are usually float32; 7 significant digits removes conversion noise.
        return float(f"{value:.7g}") if math.isfinite(value) else None
    return str(value)


def _info_value(value: object) -> object:
    """Unify INFO values: cyvcf2 returns a scalar where pysam and vcfpy return a one-item list."""

    converted = _json_value(value)
    if isinstance(converted, list) and len(converted) == 1:
        return converted[0]
    return converted


def _genotype_text(alleles: Iterable[object], phased: bool) -> str:
    text = []
    for allele in alleles:
        if allele is None or allele == "." or (isinstance(allele, int) and allele < 0):
            text.append(".")
        else:
            text.append(str(int(allele)))
    return ("|" if phased else "/").join(text)


def _site(
    *,
    chrom: object,
    pos: object,
    identifier: object,
    ref: object,
    alt: Iterable[object],
    qual: object,
    filters: Iterable[object],
    info: Mapping[str, object],
    genotypes: list[str],
) -> dict[str, Any]:
    identity = None if identifier in (None, "", ".") else str(identifier)
    return {
        "chrom": str(chrom),
        "pos": int(pos),
        "id": identity,
        "ref": str(ref),
        "alt": [str(item) for item in alt],
        "qual": _json_value(qual),
        "filter": [str(item) for item in filters if str(item) not in _IGNORED_FILTERS],
        "info": {str(key): _info_value(value) for key, value in info.items()},
        "genotypes": genotypes,
    }


def _cyvcf2_site(variant: Any) -> dict[str, Any]:
    raw_filter = variant.FILTER or ""
    return _site(
        chrom=variant.CHROM,
        pos=variant.POS,
        identifier=variant.ID,
        ref=variant.REF,
        alt=variant.ALT or [],
        qual=variant.QUAL,
        filters=str(raw_filter).split(";"),
        info=dict(variant.INFO),
        genotypes=[_genotype_text(call[:-1], bool(call[-1])) for call in variant.genotypes],
    )


def _pysam_site(variant: Any, samples: Sequence[str]) -> dict[str, Any]:
    genotypes = []
    for name in samples:
        sample = variant.samples[name]
        genotypes.append(_genotype_text(sample.get("GT") or (None,), bool(sample.phased)))
    return _site(
        chrom=variant.chrom,
        pos=variant.pos,
        identifier=variant.id,
        ref=variant.ref,
        alt=variant.alts or (),
        qual=variant.qual,
        filters=list(variant.filter.keys()),
        info=dict(variant.info.items()),
        genotypes=genotypes,
    )


def _vcfpy_site(variant: Any) -> dict[str, Any]:
    genotypes = []
    for call in variant.calls:
        alleles = call.gt_alleles if call.gt_alleles is not None else [None]
        genotypes.append(_genotype_text(alleles, bool(call.is_phased)))
    return _site(
        chrom=variant.CHROM,
        pos=variant.POS,
        identifier=";".join(variant.ID),
        ref=variant.REF,
        alt=[allele.value for allele in variant.ALT],
        qual=variant.QUAL,
        filters=variant.FILTER or [],
        info=dict(variant.INFO),
        genotypes=genotypes,
    )


def _open_reader(
    path: Path, backend: str, region: tuple[str, int | None, int | None] | None
) -> tuple[list[str], Iterator[dict[str, Any]], Any]:
    if backend == "cyvcf2":
        module = require("cyvcf2", _EXTRAS[backend])
        reader = module.VCF(str(path))
        samples = [str(name) for name in reader.samples]
        chrom, start, end = region or (None, None, None)
        if region is None:
            raw = reader
        else:
            if chrom not in reader.seqnames:
                raise ConfigurationError(
                    f"cyvcf2 could not open {path.name}: invalid contig `{chrom}`"
                )
            raw = reader(chrom if start is None else f"{chrom}:{start}-{end}")
        return samples, (_cyvcf2_site(variant) for variant in raw), reader
    if backend == "pysam":
        module = require("pysam", _EXTRAS[backend])
        handle = module.VariantFile(str(path))
        samples = [str(name) for name in handle.header.samples]
        if region is None:
            raw = handle
        else:
            chrom, start, end = region
            raw = handle.fetch(chrom) if start is None else handle.fetch(chrom, start - 1, end)
        return samples, (_pysam_site(variant, samples) for variant in raw), handle
    module = require("vcfpy", _EXTRAS[backend])
    reader = module.Reader.from_path(str(path))
    samples = [str(name) for name in reader.header.samples.names]
    return samples, (_vcfpy_site(variant) for variant in reader), reader


def _read_sites(
    path: Path,
    backend: str,
    region: tuple[str, int | None, int | None] | None,
    limit: int,
) -> tuple[list[str], list[dict[str, Any]]]:
    try:
        samples, sites, handle = _open_reader(path, backend, region)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"{backend} could not open {path.name}: {exc}") from exc
    records: list[dict[str, Any]] = []
    try:
        for site in sites:
            if len(records) >= limit:
                raise ConfigurationError(
                    f"VCF input exceeds max_records={limit}; increase the explicit limit."
                )
            records.append(site)
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"{backend} VCF parsing failed: {exc}") from exc
    finally:
        close = getattr(handle, "close", None)
        if callable(close):
            close()
    return samples, records


def _document(
    operation: str,
    path: Path,
    backend: str,
    region: str | None,
    samples: list[str],
    variants: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "backend": backend,
        "source": str(path),
        "region": region,
        "samples": samples,
        "variants": variants,
        "provenance": {
            "package": backend,
            "coordinate_system": "VCF 1-based POS",
            "genotype_convention": (
                "allele indices joined by '/' (unphased) or '|' (phased); '.' is missing"
            ),
            "filter_convention": "failing filter names only; an empty list means PASS or no filter",
        },
    }


def read_vcf(
    source: str | Path,
    *,
    backend: str = "cyvcf2",
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a VCF into dnakit variant records through one optional backend.

    ``backend`` is ``"cyvcf2"``, ``"pysam"`` or ``"vcfpy"``. The returned
    ``variants`` list can be passed directly to the native ``calculate_*``
    functions in :mod:`dnakit.variants`.
    """

    path = input_path(source)
    selected = _choose(backend, _IO_BACKENDS)
    limit = _limit(max_records)
    samples, variants = _read_sites(path, selected, None, limit)
    return _document("read_vcf", path, selected, None, samples, variants)


def fetch_vcf_region(
    source: str | Path,
    region: str,
    *,
    backend: str = "pysam",
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read only the records in ``region`` from an indexed, bgzipped VCF.

    ``region`` is 1-based and inclusive, for example ``"chr1:1000-2000"``. The
    file must have a tabix (``.tbi``) or CSI index. ``backend`` is ``"pysam"``
    or ``"cyvcf2"``.
    """

    path = input_path(source)
    selected = _choose(backend, _REGION_BACKENDS)
    parsed = _parse_region(region)
    limit = _limit(max_records)
    samples, variants = _read_sites(path, selected, parsed, limit)
    return _document("fetch_vcf_region", path, selected, region, samples, variants)


def _sample_indices(samples: Sequence[str], names: Sequence[str]) -> list[int]:
    lookup = {name: index for index, name in enumerate(samples)}
    missing = [str(name) for name in names if name not in lookup]
    if missing:
        raise ConfigurationError(f"unknown sample names: {missing}")
    return [lookup[name] for name in names]


def _names(value: object, name: str) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of sample names.")
    return [str(item) for item in value]


def _parse_call(call: object) -> list[int]:
    if not isinstance(call, str):
        raise ConfigurationError(f"genotype calls must be strings such as '0/1'; got {call!r}.")
    alleles: list[int] = []
    for token in re.split(r"[/|]", call):
        if token in ("", "."):
            alleles.append(-1)
        elif token.isdigit():
            alleles.append(int(token))
        else:
            raise ConfigurationError(f"unsupported genotype call {call!r}.")
    return alleles


def _panel(
    data: object,
    *,
    biallelic: bool,
    complete: bool,
    single_chromosome: bool,
) -> _Panel:
    """Parse variant records into a sorted scikit-allel-ready genotype array."""

    np = require("numpy", _NUMPY_EXTRA)
    records, metadata = _records(data)
    samples = _sample_names(metadata, records)
    skipped = {"multiallelic": 0, "missing_calls": 0}
    kept: list[tuple[str, int, list[list[int]]]] = []
    ploidy: int | None = None
    for record in records:
        chrom, pos, _ref, alts = _base(record)
        calls = record.get("genotypes")
        if (
            isinstance(calls, (str, bytes))
            or not isinstance(calls, Sequence)
            or len(calls) != len(samples)
        ):
            raise ConfigurationError("each variant needs one genotype call per sample.")
        parsed = [_parse_call(call) for call in calls]
        if ploidy is None:
            ploidy = len(parsed[0])
        if any(len(alleles) != ploidy for alleles in parsed):
            raise ConfigurationError("all genotype calls must share one ploidy.")
        if biallelic and len(alts) != 1:
            skipped["multiallelic"] += 1
            continue
        if complete and any(-1 in alleles for alleles in parsed):
            skipped["missing_calls"] += 1
            continue
        kept.append((chrom, pos, parsed))
    if not kept:
        raise ConfigurationError("no variant records remain after applying the site filters.")
    kept.sort(key=lambda item: (item[0], item[1]))
    chromosomes = {chrom for chrom, _pos, _calls in kept}
    if single_chromosome and len(chromosomes) != 1:
        raise ConfigurationError("this statistic needs records from exactly one chromosome.")
    return _Panel(
        samples=samples,
        chromosome=next(iter(chromosomes)) if len(chromosomes) == 1 else None,
        positions=np.array([pos for _chrom, pos, _calls in kept], dtype="int64"),
        genotypes=np.array([calls for _chrom, _pos, calls in kept], dtype="int16"),
        skipped=skipped,
    )


def _sample_names(metadata: Mapping[str, Any], records: Sequence[Mapping[str, Any]]) -> list[str]:
    if "samples" in metadata:
        return _names(metadata["samples"], "samples") if metadata["samples"] else []
    if not records:
        return []
    width = len(records[0].get("genotypes") or [])
    return [f"sample_{index + 1}" for index in range(width)]


def _selected(panel: _Panel, samples: Sequence[str] | None) -> list[int]:
    if samples is None:
        return list(range(len(panel.samples)))
    return _sample_indices(panel.samples, _names(samples, "samples"))


def _groups(panel: _Panel, populations: Mapping[str, Sequence[str]] | None) -> dict[str, list[int]]:
    if populations is None:
        return {"all": list(range(len(panel.samples)))}
    if not isinstance(populations, Mapping) or not populations:
        raise ConfigurationError(
            "populations must map each population name to a list of sample names."
        )
    return {
        str(name): _sample_indices(panel.samples, _names(members, f"population {name!r}"))
        for name, members in populations.items()
    }


def _finite(value: object) -> float | None:
    number = float(value)
    return number if math.isfinite(number) else None


def _block_length(variants: int) -> int:
    return max(1, math.ceil(variants / 10))


def _allel_genotypes(panel: _Panel, indices: Sequence[int]) -> Any:
    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    return allel.GenotypeArray(panel.genotypes).take(list(indices), axis=1)


def _require_diploid(panel: _Panel) -> None:
    if panel.genotypes.shape[2] != 2:
        raise ConfigurationError("this statistic supports diploid genotype calls only.")


def calculate_observed_heterozygosity(
    data: object,
    *,
    populations: Mapping[str, Sequence[str]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean observed heterozygosity per population, using scikit-allel.

    Each site contributes the fraction of called samples that are heterozygous.
    Sites where every call is missing are excluded from the mean.
    """

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    panel = _panel(data, biallelic=False, complete=False, single_chromosome=False)
    genotypes = allel.GenotypeArray(panel.genotypes)
    summary: dict[str, Any] = {}
    for name, indices in _groups(panel, populations).items():
        values = allel.heterozygosity_observed(genotypes.take(indices, axis=1))
        usable = [float(value) for value in values if math.isfinite(float(value))]
        summary[name] = {
            "samples": len(indices),
            "variants": len(usable),
            "mean_observed_heterozygosity": (sum(usable) / len(usable)) if usable else None,
        }
    return {
        "status": "completed",
        "operation": "calculate_observed_heterozygosity",
        "populations": summary,
        "variants_total": int(panel.positions.shape[0]),
        "skipped": panel.skipped,
        "provenance": {
            "package": "scikit-allel",
            "statistic": "per-site observed heterozygosity, averaged",
        },
    }


def calculate_fst(
    data: object,
    *,
    populations: Mapping[str, Sequence[str]],
    method: str = "hudson",
    verbose: bool = True,
) -> dict[str, Any]:
    """Genome-wide Fst between two populations, using scikit-allel.

    ``method`` is ``"hudson"`` (ratio of averages) or ``"weir_cockerham"``.
    Only biallelic sites are used. The standard error comes from a block
    jackknife with about ten blocks.
    """

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    if method not in ("hudson", "weir_cockerham"):
        raise ConfigurationError("method must be 'hudson' or 'weir_cockerham'.")
    panel = _panel(data, biallelic=True, complete=False, single_chromosome=False)
    groups = _groups(panel, populations)
    if len(groups) != 2:
        raise ConfigurationError("calculate_fst needs exactly two populations.")
    (name_a, indices_a), (name_b, indices_b) = groups.items()
    if method == "weir_cockerham" and min(len(indices_a), len(indices_b)) < 2:
        raise ConfigurationError("Weir-Cockerham Fst needs at least two samples per population.")
    genotypes = allel.GenotypeArray(panel.genotypes)
    block = _block_length(int(panel.positions.shape[0]))
    if method == "hudson":
        counts_a = genotypes.take(indices_a, axis=1).count_alleles()
        counts_b = genotypes.take(indices_b, axis=1).count_alleles()
        fst, standard_error, _blocks, _jackknife = allel.average_hudson_fst(
            counts_a, counts_b, block
        )
    else:
        fst, standard_error, _blocks, _jackknife = allel.average_weir_cockerham_fst(
            genotypes, [indices_a, indices_b], block
        )
    return {
        "status": "completed",
        "operation": "calculate_fst",
        "method": method,
        "population_a": {"name": name_a, "samples": len(indices_a)},
        "population_b": {"name": name_b, "samples": len(indices_b)},
        "variants": int(panel.positions.shape[0]),
        "block_size": block,
        "fst": _finite(fst),
        "standard_error": _finite(standard_error),
        "skipped": panel.skipped,
        "provenance": {"package": "scikit-allel", "statistic": f"{method} Fst"},
    }


def calculate_pca(
    data: object,
    *,
    n_components: int = 10,
    samples: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Principal components of sample genotypes, using scikit-allel.

    Uses biallelic, diploid sites with complete calls. Monomorphic sites are
    removed before the Patterson scaling, which is scikit-allel's default.
    """

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    require("scipy", _SCIKIT_ALLEL_EXTRA)
    panel = _panel(data, biallelic=True, complete=True, single_chromosome=False)
    _require_diploid(panel)
    indices = _selected(panel, samples)
    if isinstance(n_components, bool) or not isinstance(n_components, int) or n_components < 1:
        raise ConfigurationError("n_components must be a positive integer.")
    genotypes = _allel_genotypes(panel, indices).to_n_alt()
    alt_counts = genotypes.sum(axis=1)
    polymorphic = (alt_counts > 0) & (alt_counts < 2 * genotypes.shape[1])
    genotypes = genotypes[polymorphic]
    if genotypes.shape[0] == 0 or n_components > min(genotypes.shape):
        raise ConfigurationError(
            "n_components cannot exceed the number of samples or polymorphic sites."
        )
    coordinates, model = allel.pca(genotypes, n_components=n_components, scaler="patterson")
    names = [panel.samples[index] for index in indices]
    return {
        "status": "completed",
        "operation": "calculate_pca",
        "components": n_components,
        "samples": names,
        "coordinates": {
            name: [float(value) for value in row]
            for name, row in zip(names, coordinates, strict=True)
        },
        "explained_variance_ratio": [float(value) for value in model.explained_variance_ratio_],
        "polymorphic_variants": int(genotypes.shape[0]),
        "skipped": panel.skipped,
        "provenance": {"package": "scikit-allel", "statistic": "PCA with Patterson scaling"},
    }


def calculate_tajimas_d(
    data: object, *, samples: Sequence[str] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Tajima's D for one set of samples on one chromosome, using scikit-allel."""

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    panel = _panel(data, biallelic=True, complete=False, single_chromosome=True)
    indices = _selected(panel, samples)
    counts = _allel_genotypes(panel, indices).count_alleles()
    value = allel.tajima_d(counts, pos=panel.positions)
    return {
        "status": "completed",
        "operation": "calculate_tajimas_d",
        "chromosome": panel.chromosome,
        "samples": len(indices),
        "variants": int(panel.positions.shape[0]),
        "tajimas_d": _finite(value),
        "skipped": panel.skipped,
        "provenance": {"package": "scikit-allel", "statistic": "Tajima's D"},
    }


def calculate_watterson_theta(
    data: object,
    *,
    sequence_length: int,
    samples: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Watterson's theta per base over ``1..sequence_length``, using scikit-allel."""

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    if (
        isinstance(sequence_length, bool)
        or not isinstance(sequence_length, int)
        or sequence_length < 1
    ):
        raise ConfigurationError("sequence_length must be a positive integer.")
    panel = _panel(data, biallelic=True, complete=False, single_chromosome=True)
    if int(panel.positions.max()) > sequence_length:
        raise ConfigurationError("sequence_length must cover every variant position.")
    indices = _selected(panel, samples)
    counts = _allel_genotypes(panel, indices).count_alleles()
    value = allel.watterson_theta(panel.positions, counts, start=1, stop=sequence_length)
    return {
        "status": "completed",
        "operation": "calculate_watterson_theta",
        "chromosome": panel.chromosome,
        "sequence_length": sequence_length,
        "samples": len(indices),
        "variants": int(panel.positions.shape[0]),
        "theta_per_base": _finite(value),
        "skipped": panel.skipped,
        "provenance": {"package": "scikit-allel", "statistic": "Watterson's theta per base"},
    }


def calculate_nucleotide_diversity(
    data: object,
    *,
    sequence_length: int,
    samples: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide diversity (pi) per base over ``1..sequence_length``, using scikit-allel."""

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    if (
        isinstance(sequence_length, bool)
        or not isinstance(sequence_length, int)
        or sequence_length < 1
    ):
        raise ConfigurationError("sequence_length must be a positive integer.")
    panel = _panel(data, biallelic=True, complete=False, single_chromosome=True)
    if int(panel.positions.max()) > sequence_length:
        raise ConfigurationError("sequence_length must cover every variant position.")
    indices = _selected(panel, samples)
    counts = _allel_genotypes(panel, indices).count_alleles()
    value = allel.sequence_diversity(panel.positions, counts, start=1, stop=sequence_length)
    return {
        "status": "completed",
        "operation": "calculate_nucleotide_diversity",
        "chromosome": panel.chromosome,
        "sequence_length": sequence_length,
        "samples": len(indices),
        "variants": int(panel.positions.shape[0]),
        "pi_per_base": _finite(value),
        "skipped": panel.skipped,
        "provenance": {"package": "scikit-allel", "statistic": "nucleotide diversity per base"},
    }


def calculate_site_frequency_spectrum(
    data: object,
    *,
    samples: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unfolded site frequency spectrum from biallelic sites, using scikit-allel.

    ``spectrum[i]`` counts sites whose alternate allele appears ``i`` times among
    ``2 * samples`` diploid chromosomes.
    """

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    panel = _panel(data, biallelic=True, complete=False, single_chromosome=False)
    _require_diploid(panel)
    indices = _selected(panel, samples)
    counts = _allel_genotypes(panel, indices).count_alleles()
    chromosomes = 2 * len(indices)
    spectrum = allel.sfs(counts[:, 1], n=chromosomes)
    return {
        "status": "completed",
        "operation": "calculate_site_frequency_spectrum",
        "chromosomes": chromosomes,
        "spectrum": [int(value) for value in spectrum],
        "variants": int(panel.positions.shape[0]),
        "skipped": panel.skipped,
        "provenance": {"package": "scikit-allel", "statistic": "unfolded site frequency spectrum"},
    }


def calculate_windowed_ld(
    data: object,
    *,
    window_size: int = 10_000,
    samples: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Median r-squared linkage disequilibrium in fixed-size windows, using scikit-allel.

    Uses biallelic, diploid sites with complete calls on one chromosome.
    Window coordinates are 1-based and inclusive.
    """

    allel = require("allel", _SCIKIT_ALLEL_EXTRA)
    if isinstance(window_size, bool) or not isinstance(window_size, int) or window_size < 1:
        raise ConfigurationError("window_size must be a positive integer.")
    panel = _panel(data, biallelic=True, complete=True, single_chromosome=True)
    _require_diploid(panel)
    indices = _selected(panel, samples)
    genotypes = _allel_genotypes(panel, indices).to_n_alt()
    # Anchor windows at position 1 so boundaries are fixed rather than set by the retained sites.
    values, windows, counts = allel.windowed_r_squared(
        panel.positions, genotypes, size=window_size, start=1, stop=int(panel.positions.max())
    )
    return {
        "status": "completed",
        "operation": "calculate_windowed_ld",
        "chromosome": panel.chromosome,
        "window_size": window_size,
        "windows": [
            {
                "start": int(start),
                "end": int(stop),
                "variants": int(count),
                "median_r_squared": _finite(value),
            }
            for (start, stop), count, value in zip(windows, counts, values, strict=True)
        ],
        "variants": int(panel.positions.shape[0]),
        "skipped": panel.skipped,
        "provenance": {
            "package": "scikit-allel",
            "statistic": "windowed median r-squared (Rogers and Huff)",
        },
    }


__all__ = [
    "calculate_fst",
    "calculate_nucleotide_diversity",
    "calculate_observed_heterozygosity",
    "calculate_pca",
    "calculate_site_frequency_spectrum",
    "calculate_tajimas_d",
    "calculate_watterson_theta",
    "calculate_windowed_ld",
    "fetch_vcf_region",
    "read_vcf",
]

_register_public_verbose(globals(), tuple(__all__))
