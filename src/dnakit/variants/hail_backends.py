"""Optional Hail adapters for VCF import, MatrixTable QC, LD pruning, PCA and association.

The functions call Hail directly. Inputs are VCF files (plain or bgzip-compressed)
or MatrixTable directories written by ``import_vcf_matrix_table()``. The analyses
use biallelic sites only, because the Hail methods used here require them; the
number of skipped multi-allelic sites is returned with each result. Phenotype
files are read in Python with an explicit row limit and joined with ``hl.literal``.

Hail is imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-hail]"`` hint. The
first call in a process starts a local Spark session with a 1 GB driver.
"""

from __future__ import annotations

import csv
import importlib.metadata
import math
import os
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, positive, require
from dnakit.exceptions import ConfigurationError

_PACKAGE = "hail"
_EXTRA = "bio-hail"
_REFERENCE_GENOMES = ("GRCh37", "GRCh38")
_MODELS = ("linear", "logistic")
_MAX_LIMIT = 1_000_000
_DEFAULT_MAX_VARIANTS = 100_000
_DEFAULT_MAX_SAMPLES = 100_000
_DRIVER_MEMORY = "1g"
_MISSING_TOKENS = frozenset({"", "NA", "."})


def _hail() -> Any:
    hl = require(_PACKAGE, _EXTRA)
    log = Path(tempfile.gettempdir()) / "dnakit-hail.log"
    _start_quietly(
        lambda: hl.init(quiet=True, log=str(log), idempotent=True, driver_memory=_DRIVER_MEMORY)
    )
    return hl


def _start_quietly(start: Any) -> None:
    """Run Hail start-up while the JVM's notices on file descriptor 1 are discarded.

    ``quiet=True`` does not cover the native-library and Spark configuration notices that the
    JVM writes to the process stdout when it starts, so they are dropped for the start only.
    """
    try:
        sys.stdout.flush()
        saved = os.dup(1)
    except OSError:
        start()
        return
    try:
        with open(os.devnull, "wb") as sink:
            os.dup2(sink.fileno(), 1)
            try:
                start()
            finally:
                os.dup2(saved, 1)
    finally:
        os.close(saved)


def _limit(value: int, name: str) -> int:
    positive(value, name)
    if value > _MAX_LIMIT:
        raise ConfigurationError(f"{name} must be at most {_MAX_LIMIT}.")
    return value


def _reference(value: str) -> str:
    if value not in _REFERENCE_GENOMES:
        raise ConfigurationError(f"reference_genome must be one of {_REFERENCE_GENOMES}.")
    return value


def _boolean(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value or "\0" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL characters.")
    return value


def _fresh_path(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output path already exists: {path}")
    return path


def _number(value: Any) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _provenance(statistic: str, **notes: str) -> dict[str, str]:
    version = importlib.metadata.version(_PACKAGE)
    return {"package": f"{_PACKAGE} {version}", "statistic": statistic, **notes}


def _load(hl: Any, source: str | Path, reference_genome: str) -> Any:
    path = Path(source).expanduser().resolve()
    if path.is_dir():
        if not (path / "metadata.json.gz").is_file():
            raise ConfigurationError(f"Directory is not a Hail MatrixTable: {path}")
        return hl.read_matrix_table(str(path))
    return hl.import_vcf(str(input_path(path)), reference_genome=reference_genome)


def _biallelic(hl: Any, mt: Any, max_variants: int) -> tuple[Any, int, int]:
    """Keep biallelic sites; return the matrix, kept count and skipped multi-allelic count."""
    total = mt.count_rows()
    kept = mt.filter_rows(hl.len(mt.alleles) == 2)
    count = kept.count_rows()
    if count > max_variants:
        raise ConfigurationError(f"{count} biallelic variants exceed max_variants={max_variants}.")
    return kept, count, total - count


def _site(locus: Any, alleles: Sequence[str]) -> dict[str, Any]:
    return {
        "chrom": str(locus.contig),
        "pos": int(locus.position),
        "ref": str(alleles[0]),
        "alt": str(alleles[1]),
    }


def _site_key(site: dict[str, Any]) -> tuple[str, int, str, str]:
    return (site["chrom"], site["pos"], site["ref"], site["alt"])


def _phenotype_value(raw: str, binary: bool, sample: str) -> float | None:
    text = raw.strip()
    if text in _MISSING_TOKENS:
        return None
    try:
        value = float(text)
    except ValueError as exc:
        raise ConfigurationError(f"Phenotype for {sample!r} is not numeric: {text!r}.") from exc
    if not math.isfinite(value):
        raise ConfigurationError(f"Phenotype for {sample!r} must be finite.")
    if binary and value not in (0.0, 1.0):
        raise ConfigurationError(f"Logistic phenotype for {sample!r} must be 0, 1 or missing.")
    return value


def _read_phenotypes(
    source: str | Path, id_column: str, column: str, max_samples: int, binary: bool
) -> dict[str, float | None]:
    """Read a tab-separated phenotype table, stopping after max_samples data rows."""
    path = input_path(source)
    values: dict[str, float | None] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader, None) or []
        if id_column not in header or column not in header:
            raise ConfigurationError(
                f"{path.name} must contain the columns {id_column!r} and {column!r}."
            )
        id_index = header.index(id_column)
        value_index = header.index(column)
        for row in reader:
            if not any(cell.strip() for cell in row):
                continue
            if len(row) <= max(id_index, value_index):
                raise ConfigurationError(f"Phenotype line {reader.line_num} is too short.")
            sample = row[id_index].strip()
            if not sample or "\0" in sample:
                raise ConfigurationError(f"Phenotype line {reader.line_num} has no valid sample.")
            if sample in values:
                raise ConfigurationError(f"Duplicate phenotype sample: {sample!r}.")
            if len(values) >= max_samples:
                raise ConfigurationError(f"Phenotype file exceeds max_samples={max_samples}.")
            values[sample] = _phenotype_value(row[value_index], binary, sample)
    if not values:
        raise ConfigurationError("Phenotype file contains no samples.")
    return values


def import_vcf_matrix_table(
    source: str | Path,
    output_path: str | Path,
    *,
    reference_genome: str = "GRCh38",
    verbose: bool = True,
) -> dict[str, Any]:
    """Import a VCF with Hail and write it as a MatrixTable directory."""
    reference = _reference(reference_genome)
    target = _fresh_path(output_path)
    hl = _hail()
    mt = hl.import_vcf(str(input_path(source)), reference_genome=reference)
    n_variants, n_samples = mt.count()
    mt.write(str(target), overwrite=False)
    return {
        "status": "completed",
        "operation": "import_vcf_matrix_table",
        "path": str(target),
        "n_variants": int(n_variants),
        "n_samples": int(n_samples),
        "provenance": _provenance("VCF import and MatrixTable write", reference_genome=reference),
    }


def summarize_vcf_matrix_table(
    source: str | Path,
    *,
    reference_genome: str = "GRCh38",
    max_contigs: int = 100,
    max_samples: int = _DEFAULT_MAX_SAMPLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count variants, allele classes, contigs and samples in a VCF or MatrixTable."""
    reference = _reference(reference_genome)
    contig_limit = _limit(max_contigs, "max_contigs")
    sample_limit = _limit(max_samples, "max_samples")
    hl = _hail()
    mt = _load(hl, source, reference)
    n_variants, n_samples = mt.count()
    if n_samples > sample_limit:
        raise ConfigurationError(f"{n_samples} samples exceed max_samples={sample_limit}.")
    stats = mt.aggregate_rows(
        hl.struct(
            biallelic=hl.agg.count_where(hl.len(mt.alleles) == 2),
            snv=hl.agg.count_where(
                (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
            ),
            contigs=hl.agg.counter(mt.locus.contig),
        )
    )
    contigs = stats.contigs
    if len(contigs) > contig_limit:
        raise ConfigurationError(f"{len(contigs)} contigs exceed max_contigs={contig_limit}.")
    samples = sorted(str(sample) for sample in mt.aggregate_cols(hl.agg.collect(mt.s)))
    return {
        "status": "completed",
        "operation": "summarize_vcf_matrix_table",
        "n_variants": int(n_variants),
        "n_samples": int(n_samples),
        "biallelic_variants": int(stats.biallelic),
        "multiallelic_variants": int(n_variants) - int(stats.biallelic),
        "snv_variants": int(stats.snv),
        "contigs": [{"contig": name, "variants": int(contigs[name])} for name in sorted(contigs)],
        "samples": samples,
        "provenance": _provenance("MatrixTable row and column counts", reference_genome=reference),
    }


def calculate_variant_qc(
    source: str | Path,
    *,
    reference_genome: str = "GRCh38",
    max_variants: int = _DEFAULT_MAX_VARIANTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute per-site call rate, allele frequency and heterozygote counts."""
    reference = _reference(reference_genome)
    limit = _limit(max_variants, "max_variants")
    hl = _hail()
    mt, kept, skipped = _biallelic(hl, _load(hl, source, reference), limit)
    rows = hl.variant_qc(mt).rows().select("variant_qc").collect()
    variants = []
    for row in rows:
        qc = row.variant_qc
        variants.append(
            {
                **_site(row.locus, row.alleles),
                "n_called": int(qc.n_called),
                "n_not_called": int(qc.n_not_called),
                "call_rate": _number(qc.call_rate),
                "allele_count": int(qc.AC[1]),
                "allele_number": int(qc.AN),
                "allele_frequency": _number(qc.AF[1]),
                "n_het": int(qc.n_het),
                "n_hom_alt": int(qc.homozygote_count[1]),
            }
        )
    variants.sort(key=_site_key)
    return {
        "status": "completed",
        "operation": "calculate_variant_qc",
        "n_variants": kept,
        "skipped_multiallelic": skipped,
        "variants": variants,
        "provenance": _provenance("hl.variant_qc on biallelic sites", reference_genome=reference),
    }


def calculate_sample_qc(
    source: str | Path,
    *,
    reference_genome: str = "GRCh38",
    max_samples: int = _DEFAULT_MAX_SAMPLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute per-sample call rate, genotype class counts and transition/transversion ratios."""
    reference = _reference(reference_genome)
    sample_limit = _limit(max_samples, "max_samples")
    hl = _hail()
    mt = _load(hl, source, reference)
    n_samples = mt.count_cols()
    if n_samples > sample_limit:
        raise ConfigurationError(f"{n_samples} samples exceed max_samples={sample_limit}.")
    rows = hl.sample_qc(mt).cols().select("sample_qc").collect()
    samples = []
    for row in rows:
        qc = row.sample_qc
        samples.append(
            {
                "sample": str(row.s),
                "call_rate": _number(qc.call_rate),
                "n_called": int(qc.n_called),
                "n_not_called": int(qc.n_not_called),
                "n_hom_ref": int(qc.n_hom_ref),
                "n_het": int(qc.n_het),
                "n_hom_var": int(qc.n_hom_var),
                "r_het_hom_var": _number(qc.r_het_hom_var),
                "n_transition": int(qc.n_transition),
                "n_transversion": int(qc.n_transversion),
                "r_ti_tv": _number(qc.r_ti_tv),
            }
        )
    samples.sort(key=lambda item: item["sample"])
    return {
        "status": "completed",
        "operation": "calculate_sample_qc",
        "n_samples": n_samples,
        "samples": samples,
        "provenance": _provenance("hl.sample_qc over all sites", reference_genome=reference),
    }


def calculate_hardy_weinberg_pvalues(
    source: str | Path,
    *,
    reference_genome: str = "GRCh38",
    one_sided: bool = False,
    max_variants: int = _DEFAULT_MAX_VARIANTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test each biallelic site for Hardy-Weinberg equilibrium with an exact mid-p test."""
    reference = _reference(reference_genome)
    flag = _boolean(one_sided, "one_sided")
    limit = _limit(max_variants, "max_variants")
    hl = _hail()
    mt, kept, skipped = _biallelic(hl, _load(hl, source, reference), limit)
    mt = mt.annotate_rows(
        n_hom_ref=hl.int32(hl.agg.count_where(mt.GT.is_hom_ref())),
        n_het=hl.int32(hl.agg.count_where(mt.GT.is_het())),
        n_hom_var=hl.int32(hl.agg.count_where(mt.GT.is_hom_var())),
    )
    mt = mt.annotate_rows(
        hwe=hl.hardy_weinberg_test(mt.n_hom_ref, mt.n_het, mt.n_hom_var, one_sided=flag)
    )
    rows = mt.rows().select("n_hom_ref", "n_het", "n_hom_var", "hwe").collect()
    variants = [
        {
            **_site(row.locus, row.alleles),
            "n_hom_ref": int(row.n_hom_ref),
            "n_het": int(row.n_het),
            "n_hom_var": int(row.n_hom_var),
            "het_freq_hwe": _number(row.hwe.het_freq_hwe),
            "p_value": _number(row.hwe.p_value),
        }
        for row in rows
    ]
    variants.sort(key=_site_key)
    alternative = "excess_heterozygosity" if flag else "two_sided"
    return {
        "status": "completed",
        "operation": "calculate_hardy_weinberg_pvalues",
        "n_variants": kept,
        "skipped_multiallelic": skipped,
        "alternative": alternative,
        "variants": variants,
        "provenance": _provenance("hl.hardy_weinberg_test exact mid-p", reference_genome=reference),
    }


def prune_linkage_disequilibrium(
    source: str | Path,
    *,
    reference_genome: str = "GRCh38",
    r2: float = 0.2,
    bp_window_size: int = 1000000,
    max_variants: int = _DEFAULT_MAX_VARIANTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a maximal set of biallelic sites with pairwise r-squared below r2 in each window."""
    reference = _reference(reference_genome)
    if isinstance(r2, bool) or not isinstance(r2, (int, float)) or not 0 < r2 <= 1:
        raise ConfigurationError("r2 must be a number in (0, 1].")
    window = positive(bp_window_size, "bp_window_size")
    limit = _limit(max_variants, "max_variants")
    hl = _hail()
    mt, kept, skipped = _biallelic(hl, _load(hl, source, reference), limit)
    pruned = hl.ld_prune(mt.GT, r2=float(r2), bp_window_size=window)
    retained = sorted((_site(row.locus, row.alleles) for row in pruned.collect()), key=_site_key)
    return {
        "status": "completed",
        "operation": "prune_linkage_disequilibrium",
        "n_input_variants": kept,
        "n_retained": len(retained),
        "skipped_multiallelic": skipped,
        "retained_variants": retained,
        "provenance": _provenance(
            "hl.ld_prune on diploid genotypes, phase ignored",
            reference_genome=reference,
            r2=repr(float(r2)),
            bp_window_size=str(window),
        ),
    }


def calculate_hwe_normalized_pca(
    source: str | Path,
    *,
    reference_genome: str = "GRCh38",
    n_components: int = 10,
    max_variants: int = _DEFAULT_MAX_VARIANTS,
    max_samples: int = _DEFAULT_MAX_SAMPLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run Hardy-Weinberg-normalized PCA on biallelic genotypes and return sample scores."""
    reference = _reference(reference_genome)
    components = positive(n_components, "n_components")
    variant_limit = _limit(max_variants, "max_variants")
    sample_limit = _limit(max_samples, "max_samples")
    hl = _hail()
    mt, kept, skipped = _biallelic(hl, _load(hl, source, reference), variant_limit)
    n_samples = mt.count_cols()
    if n_samples > sample_limit:
        raise ConfigurationError(f"{n_samples} samples exceed max_samples={sample_limit}.")
    if components > n_samples:
        raise ConfigurationError("n_components cannot exceed the number of samples.")
    eigenvalues, scores, _ = hl.hwe_normalized_pca(mt.GT, k=components, compute_loadings=False)
    samples = [
        {"sample": str(row.s), "scores": [_number(value) for value in row.scores]}
        for row in scores.collect()
    ]
    samples.sort(key=lambda item: item["sample"])
    return {
        "status": "completed",
        "operation": "calculate_hwe_normalized_pca",
        "n_variants": kept,
        "n_samples": n_samples,
        "skipped_multiallelic": skipped,
        "eigenvalues": [_number(value) for value in eigenvalues],
        "samples": samples,
        "provenance": _provenance(
            "hl.hwe_normalized_pca; monomorphic sites removed by Hail",
            reference_genome=reference,
        ),
    }


def calculate_association_tests(
    source: str | Path,
    phenotype_path: str | Path,
    *,
    phenotype_column: str,
    model: str = "linear",
    sample_column: str = "sample",
    reference_genome: str = "GRCh38",
    max_samples: int = _DEFAULT_MAX_SAMPLES,
    max_variants: int = _DEFAULT_MAX_VARIANTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Intercept-only linear or logistic association test for each biallelic site."""
    reference = _reference(reference_genome)
    if model not in _MODELS:
        raise ConfigurationError(f"model must be one of {_MODELS}.")
    column = _text(phenotype_column, "phenotype_column")
    id_column = _text(sample_column, "sample_column")
    sample_limit = _limit(max_samples, "max_samples")
    variant_limit = _limit(max_variants, "max_variants")
    phenotypes = _read_phenotypes(
        phenotype_path, id_column, column, sample_limit, binary=model == "logistic"
    )
    hl = _hail()
    mt, kept, skipped = _biallelic(hl, _load(hl, source, reference), variant_limit)
    if mt.count_cols() > sample_limit:
        raise ConfigurationError(f"VCF samples exceed max_samples={sample_limit}.")
    table = hl.literal(phenotypes, dtype=hl.tdict(hl.tstr, hl.tfloat64))
    mt = mt.annotate_cols(pheno=table[mt.s])
    n_with_phenotype = int(mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.pheno))))
    if model == "linear":
        fitted = hl.linear_regression_rows(y=mt.pheno, x=mt.GT.n_alt_alleles(), covariates=[1.0])
        statistic = "t"
    else:
        fitted = hl.logistic_regression_rows(
            test="wald", y=mt.pheno, x=mt.GT.n_alt_alleles(), covariates=[1.0]
        )
        statistic = "z"
    results = []
    for row in fitted.collect():
        item: dict[str, Any] = {
            **_site(row.locus, row.alleles),
            "beta": _number(row.beta),
            "standard_error": _number(row.standard_error),
            "statistic": _number(row.t_stat if model == "linear" else row.z_stat),
            "p_value": _number(row.p_value),
        }
        if model == "logistic":
            item["converged"] = bool(row.fit.converged)
        results.append(item)
    results.sort(key=lambda item: _site_key(item))
    return {
        "status": "completed",
        "operation": "calculate_association_tests",
        "model": model,
        "phenotype_column": column,
        "covariates": [],
        "n_samples_with_phenotype": n_with_phenotype,
        "n_variants": kept,
        "skipped_multiallelic": skipped,
        "statistic": statistic,
        "results": results,
        "provenance": _provenance(
            f"hail {model} regression, intercept only",
            reference_genome=reference,
        ),
    }


__all__ = [
    "calculate_association_tests",
    "calculate_hardy_weinberg_pvalues",
    "calculate_hwe_normalized_pca",
    "calculate_sample_qc",
    "calculate_variant_qc",
    "import_vcf_matrix_table",
    "prune_linkage_disequilibrium",
    "summarize_vcf_matrix_table",
]

_register_public_verbose(globals(), tuple(__all__))
