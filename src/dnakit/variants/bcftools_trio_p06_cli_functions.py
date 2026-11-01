"""bcftools trio and family plugin facades (pack P06): Mendelian consistency, de novo screening, parental origin.

Each public function runs one fixed ``bcftools +plugin`` command (bcftools 1.24 is the tested version)
without a shell and with an explicit argument list; the user installs bcftools with its plugins (conda,
apt or a release) and puts it on PATH or passes ``bcftools_path``. A trio is given either as a PED file
(``ped_path``: family, child, father, mother, sex) or as ``trio=[child, father, mother]``. Outputs are written
to new explicit paths (an existing file is never overwritten); text reports are saved to ``output_path``
and parsed into ``report``, record-changing plugins write a VCF/BCF. Results are JSON-serialisable
dictionaries with status, backend, version, command, outputs and a bounded log tail.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics import _htslib_p06 as h
from dnakit.exceptions import ConfigurationError

_MENDELIAN_MODES = frozenset("adeEgmMS")
_DNM_MODELS = {"naive": "--use-NAIVE", "dng": "--use-DNG", "alm": "--use-ALM", "dmm": "--use-DMM"}
_CNV_TYPES = ("del", "dup")
_DNM_TAG = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?::(?:log|flag|phred|prob))?$")
_TRIO_PREFIX = re.compile(r"^(?:[12]X:)?[^,\s:]+,[^,\s:]+,[^,\s:]+$")


def _family(
    ped_path: str | Path | None,
    trio: Sequence[str] | None,
    *,
    ped_flag: str,
    trio_flag: str,
    allow_x_prefix: bool = False,
) -> list[str]:
    """Options naming the trio(s): exactly one of a PED file and ``[child, father, mother]``."""
    if (ped_path is None) == (trio is None):
        raise ConfigurationError("Give exactly one of ped_path and trio.")
    if ped_path is not None:
        return [ped_flag, str(cli.existing_file(ped_path, "ped_path"))]
    names = h.text_list(list(trio or []), "trio")
    if len(names) != 3:
        raise ConfigurationError("trio must be [child, father, mother].")
    value = ",".join(names)
    if allow_x_prefix and not _TRIO_PREFIX.match(value):
        raise ConfigurationError("trio names must not contain ':' or spaces.")
    return [trio_flag, value]


# ---------------------------------------------------------------------------
# Public operations
# ---------------------------------------------------------------------------


def bcftools_parental_origin(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    trio: Sequence[str],
    cnv_type: str,
    region: str,
    min_binomial_probability: float = 1e-2,
    greedy: bool = False,
    list_informative_sites: bool = False,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Determine the parental origin of a CNV region in a trio (``bcftools +parental-origin``).

    ``trio`` is ``[proband, father, mother]`` (``-p``), ``cnv_type`` is ``del`` or ``dup`` (``-t``) and
    ``region`` the CNV as ``chr:start-end`` (``-r``). Parental heterozygous sites with a skewed ALT
    fraction are excluded below ``min_binomial_probability`` (``-b``); ``greedy`` also uses ambiguous sites
    (``-g``); ``list_informative_sites`` lists them (``-d``). The report (type, predicted origin, quality,
    number of markers) is saved to ``output_path``; ``predicted_origin`` is the second column.
    """
    names = h.text_list(list(trio), "trio")
    if len(names) != 3:
        raise ConfigurationError("trio must be [proband, father, mother].")
    options = [
        "-p", ",".join(names),
        "-t", cli.choice(cnv_type, _CNV_TYPES, "cnv_type"),
        "-r", cli.text(region, "region"),
        "-b", repr(cli.number(min_binomial_probability, "min_binomial_probability", minimum=0.0, maximum=1.0)),
    ]
    if h.flag(greedy, "greedy"):
        options.append("-g")
    if h.flag(list_informative_sites, "list_informative_sites"):
        options.append("-d")
    if include is not None and exclude is not None:
        raise ConfigurationError("include and exclude are mutually exclusive.")
    if include is not None:
        options += ["-i", cli.text(include, "include")]
    if exclude is not None:
        options += ["-e", cli.text(exclude, "exclude")]
    run = h.run_text_plugin(
        "bcftools_parental_origin", "parental-origin", vcf_path, output_path, options=options,
        separator=False, header="none", bcftools_path=bcftools_path,
    )
    for row in run["report"]["rows"]:
        if len(row) >= 4 and row[0] in _CNV_TYPES:
            run["predicted_origin"] = row[1]
            run["quality"] = row[2]
            run["markers"] = row[3]
            break
    return run


def bcftools_trio_stats(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    ped_path: str | Path | None = None,
    trio: Sequence[str] | None = None,
    max_alt_trios: int = 0,
    debug: Sequence[str] | None = None,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transmission rate and Mendelian statistics in trio children (``bcftools +trio-stats``).

    The trios come from ``ped_path`` (``-p``) or ``trio=[proband, father, mother]`` (``-P``). The report
    has one ``FLT`` row per filter and trio: valid trio genotypes, non-reference trio GTs, DNMs or
    Mendelian errors, novel singletons, untransmitted and transmitted singletons, transitions,
    transversions, ts/tv, homozygous and recurrent DNMs; ``trio_columns`` names the columns.
    ``include`` / ``exclude`` may hold curly-bracket ranges that are scanned together, for example
    ``GQ>{10,20,30}``. ``max_alt_trios`` limits the transmission rate to sites with at most that many
    alternate trios (``-a``, 0 = unlimited); ``debug`` is a list from ``mendel-errors`` and
    ``transmitted`` (``-d``).
    """
    options = _family(ped_path, trio, ped_flag="-p", trio_flag="-P")
    options += ["-a", str(cli.integer(max_alt_trios, "max_alt_trios", minimum=0))]
    if debug is not None:
        items = h.text_list(list(debug), "debug")
        for item in items:
            cli.choice(item, ("mendel-errors", "transmitted"), "debug")
        options += ["-d", ",".join(items)]
    run = h.run_text_plugin(
        "bcftools_trio_stats", "trio-stats", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        header="none", bcftools_path=bcftools_path,
    )
    run["trio_columns"] = [
        "filter", "child", "father", "mother", "valid_gts", "non_ref_gts", "dnm_or_mendel_errors",
        "novel_singletons", "untransmitted_singletons", "transmitted_singletons", "ts", "tv",
        "ts_tv", "homozygous_dnm", "recurrent_dnm",
    ]
    run["trios"] = [
        dict(zip(run["trio_columns"], row, strict=False))
        for row in run["report"]["rows"]
        if row and re.match(r"^FLT\d+$", str(row[0]))
    ]
    return run


def bcftools_trio_switch_rate(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    ped_path: str | Path,
    regions: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phase switch rate in trio children from phased genotypes (``bcftools +trio-switch-rate``).

    ``ped_path`` lists the trios; an optional 7th column groups the results by population. The report
    has one ``TRIO`` row per trio (father, mother, child, genotypes tested, Mendelian errors, switches,
    switch percentage) and ``POP`` rows per population; they are parsed into ``trios``.
    """
    run = h.run_text_plugin(
        "bcftools_trio_switch_rate", "trio-switch-rate", vcf_path, output_path,
        options=["-p", str(cli.existing_file(ped_path, "ped_path"))],
        general=h.selection(regions=regions), header="none", bcftools_path=bcftools_path,
    )
    names = ["father", "mother", "child", "tested", "mendelian_errors", "switches", "switch_percent"]
    run["trios"] = [
        dict(zip(names, row[1:], strict=False))
        for row in run["report"]["rows"]
        if row and row[0] == "TRIO"
    ]
    return run


def bcftools_mendelian_consistency(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    ped_path: str | Path | None = None,
    trio: Sequence[str] | None = None,
    mode: str = "c",
    rules: str = "GRCh37",
    rules_file: str | Path | None = None,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count, annotate or filter genotypes by Mendelian consistency in trios (``bcftools +mendelian2``).

    The trios come from ``ped_path`` (``-P``: columns ignored, proband, father, mother, sex 1/2) or
    ``trio=[child, father, mother]`` (``-p``; prefix the child with ``1X:`` for the male chrX pattern in
    ``trio[0]``). ``mode`` (``-m``): ``c`` (default) writes a text summary with the number of good,
    erroneous and missing genotypes per trio to ``output_path`` (parsed into ``report``); otherwise it is
    a string of VCF-output letters: ``a`` add ``INFO/MERR`` and related counts, ``d`` delete genotypes
    of inconsistent trios, ``e`` keep only sites with an error, ``E`` drop sites with an error, ``g`` keep
    sites with a good trio, ``m`` keep sites with missing genotypes, ``M`` drop them, ``S`` drop skipped
    sites; the VCF is written to ``output_path``. ``rules`` is a predefined inheritance rule set
    (``GRCh37``, ``GRCh38``, ...) and ``rules_file`` a custom file.
    """
    if not isinstance(mode, str) or not mode:
        raise ConfigurationError("mode must be 'c' or a string of the letters a d e E g m M S.")
    options = _family(ped_path, trio, ped_flag="-P", trio_flag="-p", allow_x_prefix=True)
    if rules_file is not None:
        options += ["--rules-file", str(cli.existing_file(rules_file, "rules_file"))]
    else:
        options += ["--rules", cli.text(rules, "rules")]
    general = h.selection(regions=regions, include=include, exclude=exclude)
    if mode == "c":
        return h.run_text_plugin(
            "bcftools_mendelian_consistency", "mendelian2", vcf_path, output_path,
            options=[*options, "-m", "c"], general=general, separator=False, header="none",
            bcftools_path=bcftools_path, mode=mode,
        )
    if not set(mode) <= _MENDELIAN_MODES:
        raise ConfigurationError("mode must be 'c' or a string of the letters a d e E g m M S.")
    return h.run_vcf_plugin(
        "bcftools_mendelian_consistency", "mendelian2", vcf_path, output_path,
        options=[*options, "-m", mode], general=general, separator=False, output_type=output_type,
        bcftools_path=bcftools_path, mode=mode,
    )


def bcftools_trio_denovo(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    ped_path: str | Path | None = None,
    trio: Sequence[str] | None = None,
    model: str = "dmm",
    min_score: float | None = None,
    dnm_tag: str | None = None,
    chrx_regions: str = "GRCh37",
    mutation_rate: float = 1e-8,
    strictly_novel: bool = False,
    noise_prior: float | None = None,
    phi: float | None = None,
    min_vaf: float | None = None,
    max_qm: float | None = None,
    strand_bias: float | None = None,
    allele_dropout: float | None = None,
    use_parental_ad: bool = False,
    use_child_ad: bool = False,
    force_ad: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Screen variants for de novo mutations in trios and write ``FORMAT/DNM`` (``bcftools +trio-dnm3``).

    ``model`` selects the scoring model: ``dmm`` (Dirichlet-multinomial with site-noise awareness, the
    default of bcftools 1.24 and later), ``alm`` (basic allele-likelihood model, the default up to 1.23),
    ``dng`` (original DeNovoGear model) or ``naive`` (genotypes only). The trios come from ``ped_path``
    (``-P``) or ``trio=[child, father, mother]`` (``-p``). ``dnm_tag`` is the output tag and score type
    (``DNM:log`` by default, ``DNM:flag`` for ``naive`` which only supports it; ``log``, ``flag``,
    ``phred`` or ``prob``); ``min_score`` omits the annotation
    below a log score (``-m``); ``chrx_regions`` the chrX inheritance regions or preset (``-X``);
    ``mutation_rate`` the prior (``--mrate``); ``strictly_novel`` scores only novel alleles when
    inheritance is violated (``-n``). Model options: ``noise_prior`` (``--np``) and ``phi`` (both dmm
    and alm), ``min_vaf`` and ``strand_bias`` (dmm and alm; ``strand_bias`` also dng; needs
    ``FORMAT/SP``), ``max_qm`` (dmm; negative runs without ``FORMAT/QM``), ``allele_dropout`` (alm),
    ``use_parental_ad`` (``--with-pAD``) and ``use_child_ad`` (dmm, ``--with-cAD``). ``force_ad``
    calculates VAF even when the AD field count is wrong (``--force-AD``). The annotated VCF/BCF is
    written to ``output_path``.
    """
    kind = cli.choice(model, tuple(_DNM_MODELS), "model")
    options = _family(ped_path, trio, ped_flag="-P", trio_flag="-p", allow_x_prefix=True)
    options.append(_DNM_MODELS[kind])
    if dnm_tag is None:
        dnm_tag = "DNM:flag" if kind == "naive" else "DNM:log"
    if not isinstance(dnm_tag, str) or not _DNM_TAG.match(dnm_tag):
        raise ConfigurationError("dnm_tag must look like DNM or DNM:log (log, flag, phred, prob).")
    options += ["--dnm-tag", dnm_tag, "-X", cli.text(chrx_regions, "chrx_regions")]
    options += ["--mrate", repr(cli.number(mutation_rate, "mutation_rate", minimum=0.0, maximum=1.0))]
    if min_score is not None:
        options += ["-m", repr(cli.number(min_score, "min_score"))]
    if h.flag(strictly_novel, "strictly_novel"):
        options.append("-n")
    if h.flag(force_ad, "force_ad"):
        options.append("--force-AD")
    if h.flag(use_parental_ad, "use_parental_ad"):
        options.append("--with-pAD")

    def only(value: object, name: str, flag_name: str, models: tuple[str, ...]) -> None:
        if value is None:
            return
        if kind not in models:
            raise ConfigurationError(f"{name} applies to the model(s) {', '.join(models)}, not {kind}.")
        options.extend([flag_name, repr(cli.number(value, name))])

    only(noise_prior, "noise_prior", "--np", ("dmm", "alm"))
    only(phi, "phi", "--phi", ("dmm", "alm"))
    only(min_vaf, "min_vaf", "--min-vaf", ("dmm", "alm"))
    only(max_qm, "max_qm", "--max-QM", ("dmm",))
    only(strand_bias, "strand_bias", "--sb", ("dmm", "alm", "dng"))
    only(allele_dropout, "allele_dropout", "--ad", ("alm",))
    if h.flag(use_child_ad, "use_child_ad"):
        if kind != "dmm":
            raise ConfigurationError("use_child_ad applies to the dmm model.")
        options.append("--with-cAD")
    return h.run_vcf_plugin(
        "bcftools_trio_denovo", "trio-dnm3", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        output_type=output_type, bcftools_path=bcftools_path, model=kind,
    )


__all__ = [
    "bcftools_mendelian_consistency",
    "bcftools_parental_origin",
    "bcftools_trio_denovo",
    "bcftools_trio_stats",
    "bcftools_trio_switch_rate",
]

_register_public_verbose(globals(), tuple(__all__))
