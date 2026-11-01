"""GEMMA command-line facades: summary-statistic variance components, eigen decomposition, G x E, mvLMM, BSLMM and prediction.

GEMMA (https://github.com/genetics-statistics/GEMMA, tested against 0.98.5) is a program for linear mixed models in genome-wide
association studies. The user installs it (``conda install -c bioconda gemma`` or the release binary) and puts ``gemma`` on PATH
or passes ``gemma_path``. Each function runs GEMMA with a fixed argument whitelist and no shell, writes into a new output folder
(``-outdir``) and returns a JSON-serialisable dictionary: status, backend, version, command, output files, parsed tables and the
numbers from GEMMA's log (``log_values``).

Genotypes are given either as PLINK binary files (``plink_prefix``: ``.bed``, ``.bim`` and ``.fam``) or as a BIMBAM mean
genotype file (``genotype_file``, optionally gzip-compressed, with ``phenotype_file`` and optional ``annotation_file``).
``dnakit.gwas`` already runs the standard univariate association tests; this module adds the operations listed above.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="gemma",
    executable="gemma",
    hint=(
        "Install GEMMA 0.98.5 or later (for example `conda install -c bioconda gemma`, or a release binary from "
        "https://github.com/genetics-statistics/GEMMA/releases) and put gemma on PATH or pass gemma_path."
    ),
    tested_version="0.98.5",
    citation="https://doi.org/10.1038/ng.2310",
    version_arguments=("-h",),
    version_pattern=r"GEMMA\s+(\d+\.\d+\.\d+)",
)
_LOG_VALUE = re.compile(r"^##\s*(.+?)\s*=\s*(.+?)\s*$")


def _genotypes(
    plink_prefix: str | None,
    genotype_file: str | None,
    phenotype_file: str | None,
    annotation_file: str | None,
    *,
    phenotype_required: bool = True,
) -> list[str]:
    if (plink_prefix is None) == (genotype_file is None):
        raise ConfigurationError("give exactly one of plink_prefix and genotype_file.")
    arguments: list[str] = []
    if plink_prefix is not None:
        prefix = Path(plink_prefix).expanduser().resolve()
        for suffix in (".bed", ".bim", ".fam"):
            if not prefix.with_name(prefix.name + suffix).is_file():
                raise ConfigurationError(f"PLINK file is missing: {prefix}{suffix}")
        arguments += ["-bfile", str(prefix)]
    else:
        arguments += ["-g", str(_cli.existing_file(genotype_file, "genotype_file"))]  # type: ignore[arg-type]
        if phenotype_file is None and phenotype_required:
            raise ConfigurationError("phenotype_file is required with a BIMBAM genotype_file (it also selects the individuals).")
    if phenotype_file is not None:
        arguments += ["-p", str(_cli.existing_file(phenotype_file, "phenotype_file"))]
    if annotation_file is not None:
        arguments += ["-a", str(_cli.existing_file(annotation_file, "annotation_file"))]
    return arguments


def _filters(miss: float | None, maf: float | None, hwe: float | None, r2: float | None, not_snp: bool) -> list[str]:
    arguments: list[str] = []
    for flag, value, low, high in (("-miss", miss, 0.0, 1.0), ("-maf", maf, 0.0, 0.5), ("-hwe", hwe, 0.0, 1.0), ("-r2", r2, 0.0, 1.0)):
        if value is not None:
            arguments += [flag, format(_cli.number(value, flag.lstrip("-"), minimum=low, maximum=high), ".12g")]
    if _cli.flag(not_snp, "not_snp"):
        arguments.append("-notsnp")
    return arguments


def _optional_file(flag: str, value: str | None, name: str) -> list[str]:
    return [] if value is None else [flag, str(_cli.existing_file(value, name))]


def _log_values(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    values: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = _LOG_VALUE.match(line)
        if match:
            key = re.sub(r"[^a-z0-9]+", "_", match.group(1).lower()).strip("_")
            raw = match.group(2)
            try:
                values[key] = [float(item) for item in raw.split()] if " " in raw else float(raw)
            except ValueError:
                values[key] = raw
    return values


def _table(path: Path, limit: int = 20000) -> list[dict[str, Any]]:
    """Rows of a whitespace-separated GEMMA table with a header line (numbers converted, ids kept as text)."""
    if not path.is_file():
        return []
    text_columns = {"rs", "chr", "allele1", "allele0"}
    rows: list[dict[str, Any]] = []
    names: list[str] | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if not fields:
            continue
        if names is None:
            names = fields
            continue
        if len(rows) >= limit:
            break
        rows.append({name: (value if name in text_columns else _cli.convert_number(value)) for name, value in zip(names, fields, strict=False)})
    return rows


def _gemma(
    arguments: Sequence[str],
    out: Path,
    prefix: str,
    expected: Sequence[str],
    gemma_path: str | None,
    timeout_seconds: float,
    *,
    tolerate_exit_abort: bool = False,
) -> dict[str, Any]:
    """Run GEMMA into ``out`` with prefix ``prefix`` and require the listed output files."""
    files = [out / f"{prefix}.{suffix}" for suffix in expected]
    command = [*arguments, "-outdir", str(out), "-o", prefix]
    accepted = (0, 134, -6) if tolerate_exit_abort else (0,)
    record = _cli.invoke(
        _TOOL, command, cwd=out, expected=files, explicit=gemma_path, timeout_seconds=timeout_seconds, accept_codes=accepted,
    )
    if record["return_code"] != 0:
        record["warnings"] = [
            "GEMMA 0.98.5 aborts with a heap error (exit code 134 or signal 6) after writing its result files; the files listed in "
            "outputs were all written and are returned."
        ]
    record["outputs"] = [str(path) for path in sorted(out.glob(f"{prefix}.*"))]
    return record


def gemma_snp_covariance_matrix(
    output_directory: str,
    *,
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    category_file: str | None = None,
    sample: int | None = None,
    output_prefix: str = "result",
    miss: float | None = None,
    maf: float | None = None,
    hwe: float | None = None,
    r2: float | None = None,
    not_snp: bool = False,
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the SNP variance-covariance (S) matrix of a genotype panel for summary-statistic variance components (``gemma -gs``).

    The S matrix of the SNPs of each category (``category_file``, GEMMA's SNP category file; one category when omitted) is
    estimated from the individuals with a non-missing phenotype (the phenotypes only select individuals), optionally from a random
    ``sample`` of them. Outputs in the new folder ``output_directory``: ``<prefix>.S.txt`` (S and its estimated errors, returned
    as ``s_matrix``), ``<prefix>.size.txt`` (the number of SNPs per category and of individuals, ``sizes``) and
    ``<prefix>.snps.txt.gz``. The SNP filters ``miss``, ``maf``, ``hwe`` and ``r2`` keep GEMMA's defaults when omitted.
    GEMMA 0.98.5 aborts with a heap error after it has written these files; the facade accepts that and says so in ``warnings``.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = _genotypes(plink_prefix, genotype_file, phenotype_file, None)
    arguments += _optional_file("-cat", category_file, "category_file")
    if sample is not None:
        arguments += ["-sample", str(_cli.integer(sample, "sample", minimum=2))]
    arguments += _filters(miss, maf, hwe, r2, not_snp) + ["-gs"]
    record = _gemma(arguments, out, output_prefix, ["S.txt", "size.txt"], gemma_path, timeout_seconds, tolerate_exit_abort=True)
    s_rows = [line.split("\t") for line in (out / f"{output_prefix}.S.txt").read_text().splitlines() if line]
    sizes = [float(x) for x in (out / f"{output_prefix}.size.txt").read_text().split()]
    return _cli.finish(
        record, "snp_variance_covariance_matrix", output_directory=str(out), s_matrix=[[float(x) for x in row] for row in s_rows],
        sizes=sizes, log_values=_log_values(out / f"{output_prefix}.log.txt"),
    )


def gemma_summary_statistic_q_vector(
    output_directory: str,
    beta_file: str,
    phenotype_file: str,
    *,
    category_file: str | None = None,
    output_prefix: str = "result",
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the q vector of a GWAS summary-statistic file for variance components (``gemma -gq``).

    ``beta_file`` is GEMMA's beta/z file: a header line, then the columns SNP id, sample size N, marginal z-score and the two
    alleles (the output ``*.assoc.txt`` of a GEMMA linear-model analysis is also accepted). ``category_file`` assigns the SNPs
    to categories. GEMMA 0.98.5 additionally needs a phenotype file to start (``phenotype_file``; its values are not used). The
    q vector is returned as ``q`` with its standard errors ``q_errors`` and the sizes; files are in the new folder
    ``output_directory``.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = ["-beta", str(_cli.existing_file(beta_file, "beta_file")), "-p", str(_cli.existing_file(phenotype_file, "phenotype_file"))]
    arguments += _optional_file("-cat", category_file, "category_file") + ["-gq"]
    record = _gemma(arguments, out, output_prefix, ["q.txt", "Vq.txt", "size.txt"], gemma_path, timeout_seconds)
    return _cli.finish(
        record, "summary_statistic_q_vector", output_directory=str(out),
        q=[float(x) for x in (out / f"{output_prefix}.q.txt").read_text().split()],
        q_errors=[float(x) for x in (out / f"{output_prefix}.Vq.txt").read_text().split()],
        sizes=[float(x) for x in (out / f"{output_prefix}.size.txt").read_text().split()],
    )


def gemma_eigen_decomposition(
    output_directory: str,
    kinship_file: str,
    *,
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    covariates_file: str | None = None,
    output_prefix: str = "result",
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Eigen-decompose a kinship (relatedness) matrix (``gemma -eigen``).

    GEMMA takes the rows and columns of ``kinship_file`` (a matrix from ``gemma -gk``) that belong to the analysed individuals,
    which it finds in the genotype and phenotype files, centres the matrix and decomposes it. The eigenvalues
    (``eigen_values``, descending) and the eigenvector matrix (``<prefix>.eigenU.txt``) are written to the new folder
    ``output_directory``; the eigenvalues are also returned.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = _genotypes(plink_prefix, genotype_file, phenotype_file, None, phenotype_required=True)
    arguments += ["-k", str(_cli.existing_file(kinship_file, "kinship_file"))] + _optional_file("-c", covariates_file, "covariates_file") + ["-eigen"]
    record = _gemma(arguments, out, output_prefix, ["eigenD.txt", "eigenU.txt"], gemma_path, timeout_seconds)
    values = [float(x) for x in (out / f"{output_prefix}.eigenD.txt").read_text().split()]
    return _cli.finish(
        record, "kinship_eigen_decomposition", output_directory=str(out), eigen_values=values, individual_count=len(values),
        log_values=_log_values(out / f"{output_prefix}.log.txt"),
    )


def gemma_summary_variance_components(
    output_directory: str,
    beta_file: str,
    *,
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    category_file: str | None = None,
    ld_score_file: str | None = None,
    method: str = "mqs_hew",
    covariates_file: str | None = None,
    sample: int | None = None,
    output_prefix: str = "result",
    miss: float | None = None,
    maf: float | None = None,
    hwe: float | None = None,
    r2: float | None = None,
    not_snp: bool = False,
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate variance components and heritability from GWAS summary statistics and a reference LD panel (``gemma -vc``).

    ``beta_file`` holds the marginal z-scores (see ``gemma_summary_statistic_q_vector``) and the genotypes (PLINK or BIMBAM) are
    those of the study or of a reference panel; the phenotype values are only used to select individuals. ``method`` is
    ``'mqs_hew'`` (MQS with HE weights, ``-vc 1``, the default) or ``'mqs_ldw'`` (MQS with LD-score weights, ``-vc 2``, which
    needs ``ld_score_file`` passed as ``-wcat``). ``category_file`` gives several variance components; ``covariates_file``
    centres the genotypes per group (for example reference sub-populations); ``sample`` uses a random subset of individuals.
    Result keys: ``log_values`` (GEMMA's ``pve estimates``, ``se(pve)``, ``sigma2 estimates``, ``enrichment`` ...), and the files
    in the new folder ``output_directory`` (``<prefix>.log.txt``, ``S.txt``, ``q.txt``, ``Vq.txt``, ``size.txt``).
    """
    out = _cli.new_directory(output_directory, "output_directory")
    chosen = _cli.choice(method, ("mqs_hew", "mqs_ldw"), "method")
    arguments = ["-beta", str(_cli.existing_file(beta_file, "beta_file"))]
    arguments += _genotypes(plink_prefix, genotype_file, phenotype_file, None)
    arguments += _optional_file("-cat", category_file, "category_file") + _optional_file("-c", covariates_file, "covariates_file")
    if chosen == "mqs_ldw":
        if ld_score_file is None:
            raise ConfigurationError("method='mqs_ldw' needs ld_score_file.")
        arguments += ["-wcat", str(_cli.existing_file(ld_score_file, "ld_score_file"))]
    elif ld_score_file is not None:
        raise ConfigurationError("ld_score_file is only used with method='mqs_ldw'.")
    if sample is not None:
        arguments += ["-sample", str(_cli.integer(sample, "sample", minimum=2))]
    arguments += _filters(miss, maf, hwe, r2, not_snp) + ["-vc", "2" if chosen == "mqs_ldw" else "1"]
    record = _gemma(arguments, out, output_prefix, ["log.txt"], gemma_path, timeout_seconds)
    return _cli.finish(record, "summary_statistic_variance_components", output_directory=str(out), method=chosen, log_values=_log_values(out / f"{output_prefix}.log.txt"))


def gemma_summary_variance_component_errors(
    output_directory: str,
    beta_file: str,
    reference_prefix: str,
    pve: float,
    *,
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    category_file: str | None = None,
    ld_score_file: str | None = None,
    method: str = "mqs_hew",
    output_prefix: str = "result",
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Standard errors of summary-statistic variance components with the asymptotic method and a reference panel (``gemma -ci``).

    ``reference_prefix`` is the output prefix, including its folder, of the earlier ``gemma_summary_variance_components`` (or
    ``gemma_snp_covariance_matrix`` plus ``gemma_summary_statistic_q_vector``) run on the reference genotypes (its ``S`` and ``q``
    files), and ``pve`` the PVE estimate of that run. The genotypes and phenotypes passed here must be those of all individuals of
    the study, because the asymptotic method computes extra summary statistics from them. ``method`` is ``'mqs_hew'`` or
    ``'mqs_ldw'`` (with ``ld_score_file``). The log values (``pve estimates``, ``se(pve)``, ``sigma2 per snp``, ``enrichment``
    and their errors) are returned in ``log_values``.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    chosen = _cli.choice(method, ("mqs_hew", "mqs_ldw"), "method")
    arguments = ["-beta", str(_cli.existing_file(beta_file, "beta_file"))]
    arguments += _genotypes(plink_prefix, genotype_file, phenotype_file, None)
    arguments += ["-ref", str(Path(reference_prefix).expanduser().resolve()), "-pve", format(_cli.number(pve, "pve"), ".12g")]
    arguments += _optional_file("-cat", category_file, "category_file")
    if chosen == "mqs_ldw":
        if ld_score_file is None:
            raise ConfigurationError("method='mqs_ldw' needs ld_score_file.")
        arguments += ["-wcat", str(_cli.existing_file(ld_score_file, "ld_score_file"))]
    arguments += ["-ci", "2" if chosen == "mqs_ldw" else "1"]
    record = _gemma(arguments, out, output_prefix, ["log.txt"], gemma_path, timeout_seconds)
    return _cli.finish(record, "summary_statistic_variance_component_errors", output_directory=str(out), method=chosen, log_values=_log_values(out / f"{output_prefix}.log.txt"))


def _lmm_arguments(
    kinship_file: str | None, eigen_values_file: str | None, eigen_vectors_file: str | None, covariates_file: str | None, test: str,
) -> list[str]:
    tests = {"wald": "1", "likelihood_ratio": "2", "score": "3", "all": "4"}
    arguments: list[str] = []
    if kinship_file is not None:
        arguments += ["-k", str(_cli.existing_file(kinship_file, "kinship_file"))]
    elif eigen_values_file is not None and eigen_vectors_file is not None:
        arguments += ["-d", str(_cli.existing_file(eigen_values_file, "eigen_values_file")), "-u", str(_cli.existing_file(eigen_vectors_file, "eigen_vectors_file"))]
    else:
        raise ConfigurationError("give kinship_file, or eigen_values_file together with eigen_vectors_file.")
    arguments += _optional_file("-c", covariates_file, "covariates_file")
    return arguments + ["-lmm", tests[_cli.choice(test, tuple(tests), "test")]]


def gemma_gxe_association(
    output_directory: str,
    plink_prefix: str,
    environment_file: str,
    *,
    kinship_file: str | None = None,
    eigen_values_file: str | None = None,
    eigen_vectors_file: str | None = None,
    phenotype_file: str | None = None,
    covariates_file: str | None = None,
    test: str = "wald",
    output_prefix: str = "result",
    miss: float | None = None,
    maf: float | None = None,
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test SNP x environment (G x E) interaction with a linear mixed model (``gemma -gxe -lmm``).

    For every SNP GEMMA fits a mixed model that controls the SNP main effect and the environmental main effect and tests the
    interaction. ``environment_file`` has one column with the environmental value of each individual; the relatedness comes from
    ``kinship_file`` (or the eigen files of ``gemma_eigen_decomposition``). ``test`` is ``'wald'``, ``'likelihood_ratio'``,
    ``'score'`` or ``'all'``. GEMMA 0.98.5 reads the genotypes of this analysis incorrectly from BIMBAM files, so PLINK binary
    files (``plink_prefix``) are required. The association table (``<prefix>.assoc.txt``: chromosome, SNP, position, allele
    frequency, interaction beta, se, p-values) is returned as ``associations``.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = _genotypes(plink_prefix, None, phenotype_file, None)
    arguments += ["-gxe", str(_cli.existing_file(environment_file, "environment_file"))] + _filters(miss, maf, None, None, False)
    arguments += _lmm_arguments(kinship_file, eigen_values_file, eigen_vectors_file, covariates_file, test)
    record = _gemma(arguments, out, output_prefix, ["assoc.txt"], gemma_path, timeout_seconds)
    rows = _table(out / f"{output_prefix}.assoc.txt")
    return _cli.finish(record, "gxe_linear_mixed_model", output_directory=str(out), association_count=len(rows), associations=rows[:1000], log_values=_log_values(out / f"{output_prefix}.log.txt"))


def gemma_multivariate_association(
    output_directory: str,
    phenotype_columns: Sequence[int],
    *,
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    annotation_file: str | None = None,
    kinship_file: str | None = None,
    eigen_values_file: str | None = None,
    eigen_vectors_file: str | None = None,
    covariates_file: str | None = None,
    output_prefix: str = "result",
    miss: float | None = None,
    maf: float | None = None,
    hwe: float | None = None,
    r2: float | None = None,
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Association of SNPs with several phenotypes at once with a multivariate linear mixed model (``gemma -lmm -n 1 2 ...``).

    ``phenotype_columns`` lists the 1-based columns of the phenotype file (or ``.fam`` file) to analyse jointly; at least two
    are needed, and GEMMA recommends fewer than ten. The relatedness comes from ``kinship_file`` (or the eigen files). GEMMA
    tests each SNP against all phenotypes with a Wald test; the table ``<prefix>.assoc.txt`` has one effect column per phenotype
    (``beta_1`` ...), their covariance matrix and one p-value, returned as ``associations``; the genetic and residual
    covariance of the null model are in ``log_values``.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    if isinstance(phenotype_columns, str) or len(phenotype_columns) < 2:
        raise ConfigurationError("phenotype_columns must list at least two phenotype columns.")
    columns = [str(_cli.integer(item, "phenotype_columns", minimum=1)) for item in phenotype_columns]
    arguments = _genotypes(plink_prefix, genotype_file, phenotype_file, annotation_file)
    arguments += _filters(miss, maf, hwe, r2, False)
    arguments += _lmm_arguments(kinship_file, eigen_values_file, eigen_vectors_file, covariates_file, "wald") + ["-n", *columns]
    record = _gemma(arguments, out, output_prefix, ["assoc.txt"], gemma_path, timeout_seconds)
    rows = _table(out / f"{output_prefix}.assoc.txt")
    return _cli.finish(
        record, "multivariate_linear_mixed_model", output_directory=str(out), phenotype_columns=[int(c) for c in columns],
        association_count=len(rows), associations=rows[:1000], log_values=_log_values(out / f"{output_prefix}.log.txt"),
    )


def gemma_bslmm(
    output_directory: str,
    *,
    model: str = "linear",
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    annotation_file: str | None = None,
    phenotype_column: int = 1,
    kinship_file: str | None = None,
    burn_in: int = 100_000,
    sampling_steps: int = 1_000_000,
    max_snps_in_model: int | None = None,
    seed: int | None = None,
    output_prefix: str = "result",
    miss: float | None = None,
    maf: float | None = None,
    hwe: float | None = None,
    r2: float | None = None,
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a Bayesian sparse linear mixed model for multi-locus association and polygenic prediction (``gemma -bslmm``).

    ``model`` is ``'linear'`` (standard BSLMM fitted by MCMC, ``-bslmm 1``), ``'ridge'`` (ridge regression / GBLUP, no MCMC,
    ``-bslmm 2``) or ``'probit'`` (probit BSLMM for 0/1 phenotypes, ``-bslmm 3``). ``phenotype_column`` selects the phenotype
    (``-n``); individuals with a missing value are left out and can be predicted with ``gemma_predict_phenotypes``.
    ``burn_in`` (``-w``) and ``sampling_steps`` (``-s``) are GEMMA's MCMC lengths (its defaults are 100000 and 1000000; use
    smaller values for small data), ``max_snps_in_model`` is ``-smax`` and ``seed`` ``-seed``. Outputs in the new folder
    ``output_directory``: ``<prefix>.param.txt`` (posterior mean alpha, beta and gamma per SNP, returned as ``snp_effects``),
    ``<prefix>.hyp.txt`` (hyper-parameter samples, summarised in ``hyper_parameter_means``: h, pve, rho, pge, pi, n_gamma),
    ``<prefix>.bv.txt``, ``<prefix>.gamma.txt`` (MCMC models only) and the log.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    number = {"linear": 1, "ridge": 2, "probit": 3}[_cli.choice(model, ("linear", "ridge", "probit"), "model")]
    arguments = _genotypes(plink_prefix, genotype_file, phenotype_file, annotation_file, phenotype_required=True)
    arguments += _filters(miss, maf, hwe, r2, False) + _optional_file("-k", kinship_file, "kinship_file")
    arguments += ["-n", str(_cli.integer(phenotype_column, "phenotype_column", minimum=1)), "-bslmm", str(number)]
    if number != 2:
        arguments += ["-w", str(_cli.integer(burn_in, "burn_in", minimum=0)), "-s", str(_cli.integer(sampling_steps, "sampling_steps", minimum=1))]
    if max_snps_in_model is not None:
        arguments += ["-smax", str(_cli.integer(max_snps_in_model, "max_snps_in_model", minimum=1))]
    if seed is not None:
        arguments += ["-seed", str(_cli.integer(seed, "seed", minimum=0))]
    expected = ["param.txt", "bv.txt", "log.txt"] + (["hyp.txt", "gamma.txt"] if number != 2 else [])
    record = _gemma(arguments, out, output_prefix, expected, gemma_path, timeout_seconds)
    effects = _table(out / f"{output_prefix}.param.txt")
    hyper = _table(out / f"{output_prefix}.hyp.txt")
    means = {key: sum(float(row[key]) for row in hyper) / len(hyper) for key in hyper[0]} if hyper else {}
    return _cli.finish(
        record, "bayesian_sparse_linear_mixed_model", output_directory=str(out), model=model, snp_effects=effects[:1000],
        snp_count=len(effects), hyper_parameter_means=means, log_values=_log_values(out / f"{output_prefix}.log.txt"),
    )


def gemma_predict_phenotypes(
    output_directory: str,
    *,
    source: str = "bslmm",
    plink_prefix: str | None = None,
    genotype_file: str | None = None,
    phenotype_file: str | None = None,
    parameter_file: str | None = None,
    log_file: str | None = None,
    breeding_value_file: str | None = None,
    kinship_file: str | None = None,
    phenotype_columns: Sequence[int] | None = None,
    probability_scale: bool = False,
    output_prefix: str = "result",
    gemma_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the missing phenotypes of individuals, from BSLMM estimates or from the kinship matrix (``gemma -predict``).

    ``source='bslmm'`` uses the output of ``gemma_bslmm`` run on the same genotype and phenotype files: ``parameter_file``
    (``<prefix>.param.txt``) and ``log_file`` (``<prefix>.log.txt``, the estimated mean), optionally ``breeding_value_file``
    (``<prefix>.bv.txt``) with ``kinship_file`` for the more formal predictor; ``probability_scale=True`` maps the values
    through the normal CDF (``-predict 2``, only for a probit BSLMM). ``source='kinship'`` predicts the missing values of the
    phenotype columns ``phenotype_columns`` (several give the BLUP across phenotypes, ``-n``) from ``kinship_file`` alone.
    Individuals with a missing phenotype get a predicted value, the others ``NA``; the values are returned as ``predictions``
    (one number or ``None`` per individual, in file order) next to ``<prefix>.prdt.txt`` in the new folder ``output_directory``.
    """
    out = _cli.new_directory(output_directory, "output_directory")
    mode = _cli.choice(source, ("bslmm", "kinship"), "source")
    arguments = _genotypes(plink_prefix, genotype_file, phenotype_file, None, phenotype_required=True)
    if mode == "bslmm":
        if parameter_file is None or log_file is None:
            raise ConfigurationError("source='bslmm' needs parameter_file and log_file.")
        arguments += ["-epm", str(_cli.existing_file(parameter_file, "parameter_file")), "-emu", str(_cli.existing_file(log_file, "log_file"))]
        arguments += _optional_file("-ebv", breeding_value_file, "breeding_value_file") + _optional_file("-k", kinship_file, "kinship_file")
        if phenotype_columns is not None:
            raise ConfigurationError("phenotype_columns applies to source='kinship'.")
    else:
        if kinship_file is None:
            raise ConfigurationError("source='kinship' needs kinship_file.")
        if parameter_file is not None or log_file is not None or breeding_value_file is not None:
            raise ConfigurationError("parameter_file, log_file and breeding_value_file apply to source='bslmm'.")
        arguments += ["-k", str(_cli.existing_file(kinship_file, "kinship_file"))]
        if phenotype_columns is not None:
            if isinstance(phenotype_columns, str) or not phenotype_columns:
                raise ConfigurationError("phenotype_columns must be a non-empty list of columns.")
            arguments += ["-n", *[str(_cli.integer(item, "phenotype_columns", minimum=1)) for item in phenotype_columns]]
    arguments += ["-predict", "2" if _cli.flag(probability_scale, "probability_scale") else "1"]
    record = _gemma(arguments, out, output_prefix, ["prdt.txt"], gemma_path, timeout_seconds)
    values: list[Any] = []
    for line in (out / f"{output_prefix}.prdt.txt").read_text().splitlines():
        fields = [None if item == "NA" else float(item) for item in line.split()]
        values.append(fields[0] if len(fields) == 1 else (fields or None))
    predicted = sum(1 for value in values if value is not None and value != [None] * (len(value) if isinstance(value, list) else 0))
    return _cli.finish(record, "predict_phenotypes", output_directory=str(out), source=mode, predictions=values, predicted_count=predicted)


__all__ = [
    "gemma_bslmm",
    "gemma_eigen_decomposition",
    "gemma_gxe_association",
    "gemma_multivariate_association",
    "gemma_predict_phenotypes",
    "gemma_snp_covariance_matrix",
    "gemma_summary_statistic_q_vector",
    "gemma_summary_variance_component_errors",
    "gemma_summary_variance_components",
]
_register_public_verbose(globals())
