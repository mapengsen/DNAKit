"""bcftools plugin facades (pack P06): per-site and per-sample VCF/BCF annotation, statistics and filters.

Each public function runs one fixed ``bcftools +plugin`` command (bcftools 1.24 is the tested version)
without a shell and with an explicit argument list; the user installs bcftools with its plugins (conda,
apt or a release) and puts it on PATH or passes ``bcftools_path``. Outputs are written to new explicit paths
or new folders (an existing path is never overwritten). Plugins that print a report have it saved to
``output_path`` and parsed into ``report`` (comment lines, column names and numeric rows); plugins that
change records write a VCF/BCF and the result gives the record count (plain VCF output only). The trio
and family plugins are in ``bcftools_trio_p06_cli_functions``; the core ``convert``, ``index`` and ``norm``
commands are in ``bcftools_p06_cli_functions``.
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics import _htslib_p06 as h
from dnakit.exceptions import ConfigurationError

_VARIANT_TYPES = ("snp", "indel")
_SET_TARGET = re.compile(r"^(?:\./\.|\./x|\.|a|q|b(?::[A-Za-z]+<[0-9.eE+-]+)?|r:[0-9.]+)$")
_SET_NEW = re.compile(r"^(?:[.0imMXpu]{1,2}|c:[0-9./|mMX.]+)$")
_PRUNE_STATS = ("r2", "LD", "RD", "count")
_NSITES_MODES = ("maxAF", "1st", "rand")
_DIRECTIONS = ("nearest", "fwd", "rev", "both")
_FIXREF_MODES = ("stats", "flip", "flip-all", "id", "ref-alt", "swap", "top")
_REMOVE_OVERLAP_OUTPUT = ("v", "z", "b", "u", "t", "tz")
_BIN_LIST = re.compile(r"^[0-9.eE+-]+(?:,[0-9.eE+-]+)*$")
_TAG2TAG = tuple(
    sorted(
        {f"{a}-to-{b}" for a in ("GL", "PL", "GP") for b in ("GL", "PL", "GP", "GT") if a != b}
        | {"QR-QA-to-QS", "LXX-to-XX", "LPL-to-PL", "LAD-to-AD"}
    )
)
_ANNOTATIONS = ("PASSOC", "FASSOC", "NASSOC", "NOVELAL", "NOVELGT", "BCSQ")
_CSQ_TAG = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*$")


def _number_text(value: object, name: str, *, minimum: float | None = None, maximum: float | None = None) -> str:
    return repr(cli.number(value, name, minimum=minimum, maximum=maximum))


def _bin_list(value: Sequence[float] | None, name: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or len(value) < 2:
        raise ConfigurationError(f"{name} must be a list of at least two numbers.")
    return [",".join(_number_text(item, name) for item in value)]


# ---------------------------------------------------------------------------
# Allele and genotype statistics (report plugins)
# ---------------------------------------------------------------------------


def bcftools_ad_bias(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    sample_pairs: Sequence[Sequence[str]],
    threshold: float = 1e-3,
    min_alt_depth: int = 1,
    min_depth: int = 0,
    variant_type: str | None = None,
    extra_columns: str | None = None,
    clean_vcf: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fisher test for a varying ALT allele fraction (``FORMAT/AD``) between sample pairs (``bcftools +ad-bias``).

    ``sample_pairs`` lists ``[sample, control]`` pairs. Sites with a p-value below ``threshold``
    (``-t``) are reported, one ``FT`` row per hit with the allele depths of both samples. ``min_alt_depth``
    / ``min_depth`` are ``-a`` / ``-d``; ``variant_type`` is ``snp`` or ``indel`` (``-v``);
    ``extra_columns`` is a ``bcftools query`` style format appended to every row (``-f``). With
    ``clean_vcf`` True (``-c``) the output is instead a VCF without the sites and ALT alleles that fail
    the threshold.
    """
    if isinstance(sample_pairs, (str, bytes)) or not sample_pairs:
        raise ConfigurationError("sample_pairs must be a non-empty list of [sample, control] pairs.")
    lines = []
    for pair in sample_pairs:
        if isinstance(pair, (str, bytes)) or len(pair) != 2:
            raise ConfigurationError("each sample pair must have exactly two sample names.")
        lines.append("\t".join(cli.text(name, "sample name") for name in pair))
    options = [
        "-t", _number_text(threshold, "threshold", minimum=0.0, maximum=1.0),
        "-a", str(cli.integer(min_alt_depth, "min_alt_depth", minimum=0)),
        "-d", str(cli.integer(min_depth, "min_depth", minimum=0)),
    ]
    if variant_type is not None:
        options += ["-v", cli.choice(variant_type, _VARIANT_TYPES, "variant_type")]
    if extra_columns is not None:
        options += ["-f", cli.text(extra_columns, "extra_columns")]
    general = h.selection(regions=regions, include=include, exclude=exclude)
    with tempfile.TemporaryDirectory(prefix="dnakit-adbias-") as folder:
        pairs_file = Path(folder) / "pairs.txt"
        pairs_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        options += ["-s", str(pairs_file)]
        if clean_vcf:
            return h.run_vcf_plugin(
                "bcftools_ad_bias", "ad-bias", vcf_path, output_path,
                options=[*options, "-c"], general=general, bcftools_path=bcftools_path,
            )
        return h.run_text_plugin(
            "bcftools_ad_bias", "ad-bias", vcf_path, output_path,
            options=options, general=general, bcftools_path=bcftools_path,
        )


def bcftools_af_distribution(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    deviation_bins: Sequence[float] | None = None,
    probability_bins: Sequence[float] | None = None,
    per_sample: bool = False,
    af_tag: str = "AF",
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """AF-deviation and genotype-probability distributions under Hardy-Weinberg (``bcftools +af-dist``).

    Only non-reference genotypes at sites with a known allele frequency are used (probabilities
    ``2*AF*(1-AF)`` for heterozygous and ``AF**2`` for homozygous alternate). ``deviation_bins`` /
    ``probability_bins`` are the bin edges (``-d`` / ``-p``; the bcftools default is 0, 0.1, ..., 1);
    ``af_tag`` is the INFO tag that holds the frequency (``-t``); ``per_sample`` adds the geometric mean
    probability per sample (``-s``). The report (``PROB_DIST``, ``DEV_DIST`` and sample rows) is saved to
    ``output_path``.
    """
    options: list[str] = []
    for flag_name, values, label in (
        ("-d", deviation_bins, "deviation_bins"),
        ("-p", probability_bins, "probability_bins"),
    ):
        for item in _bin_list(values, label):
            options += [flag_name, item]
    options += ["-t", cli.text(af_tag, "af_tag")]
    if h.flag(per_sample, "per_sample"):
        options.append("-s")
    return h.run_text_plugin(
        "bcftools_af_distribution", "af-dist", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude),
        header="none", bcftools_path=bcftools_path,
    )


def bcftools_allele_length_counts(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Frequency of REF and ALT allele lengths (``bcftools +allele-length``).

    The table has the columns ``LENGTH``, ``REF``, ``ALT``, ``REF+ALT`` and ``REF+ALT WITH NON-BASE
    NUCLEOTIDES``; it is saved to ``output_path`` and parsed into ``report``.
    """
    return h.run_text_plugin(
        "bcftools_allele_length_counts", "allele-length", vcf_path, output_path,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        header="first", bcftools_path=bcftools_path,
    )


def bcftools_check_ploidy(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    use_missing: bool = False,
    regions: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report the ploidy of every sample along each chromosome (``bcftools +check-ploidy``).

    One row per sample and region: sample, chromosome, region start, region end, ploidy. Missing and
    half-missing genotypes are ignored unless ``use_missing`` (``-m``). The report is saved to
    ``output_path``.
    """
    options = ["-m"] if h.flag(use_missing, "use_missing") else []
    return h.run_text_plugin(
        "bcftools_check_ploidy", "check-ploidy", vcf_path, output_path, options=options,
        general=h.selection(regions=regions), bcftools_path=bcftools_path,
    )


def bcftools_check_sparsity(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    min_markers: int = 1,
    regions: str | None = None,
    regions_file: str | Path | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the samples that have fewer than ``min_markers`` genotypes in a region or chromosome
    (``bcftools +check-sparsity``).

    Without ``regions`` / ``regions_file`` every chromosome is checked (``-n``, ``-r``, ``-R``). The
    sample list is saved to ``output_path``; an empty file means every sample has enough markers.
    """
    options = ["-n", str(cli.integer(min_markers, "min_markers", minimum=1))]
    if regions is not None:
        options += ["-r", cli.text(regions, "regions")]
    if regions_file is not None:
        options += ["-R", str(cli.existing_file(regions_file, "regions_file"))]
    return h.run_text_plugin(
        "bcftools_check_sparsity", "check-sparsity", vcf_path, output_path,
        options=options, separator=False, header="none", bcftools_path=bcftools_path,
    )


def bcftools_genotype_dosage(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    tags: Sequence[str] = ("PL", "GL", "GT"),
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype dosage per site and sample from PL, GL or GT (``bcftools +dosage``).

    ``tags`` is the priority list of FORMAT tags the dosage is taken from (``-t``, default PL, GL, GT).
    One row per site: CHROM, POS, REF, ALT, then one dosage (0 to 2, -1 for missing) per sample. The table
    is saved to ``output_path``.
    """
    options = ["-t", ",".join(h.text_list(list(tags), "tags"))]
    return h.run_text_plugin(
        "bcftools_genotype_dosage", "dosage", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude),
        bcftools_path=bcftools_path,
    )


def bcftools_genotype_intersections(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    include_missing: bool = False,
    annotate_subsets: bool = False,
    human_readable: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count genotype intersections across all subsets of samples (``bcftools +GTisec``).

    ``include_missing`` adds the missing-genotype count per sample (``-m``), ``annotate_subsets`` names
    the samples of every subset in each count row (``-v``) and ``human_readable`` sorts the output by
    sample (``-H``, implies ``annotate_subsets``). The counts are saved to ``output_path``.
    """
    options = [
        option
        for enabled, option, name in (
            (include_missing, "-m", "include_missing"),
            (annotate_subsets, "-v", "annotate_subsets"),
            (human_readable, "-H", "human_readable"),
        )
        if h.flag(enabled, name)
    ]
    return h.run_text_plugin(
        "bcftools_genotype_intersections", "GTisec", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude), header="none",
        bcftools_path=bcftools_path,
    )


def bcftools_indel_stats(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    ped_path: str | Path | None = None,
    csq_tag: str = "CSQ",
    max_length: int = 20,
    alt_to_ref_dnm: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-sample indel statistics, or de novo indel statistics for the trios of a PED file
    (``bcftools +indel-stats``).

    ``max_length`` is the longest indel counted (``--max-len``), ``csq_tag`` the INFO tag with VEP/BCSQ consequences used to separate
    frameshift and inframe indels (``-c``). With ``ped_path`` (``-p``) only de novo indels are counted;
    ``alt_to_ref_dnm`` also counts GT errors such as 0/1 + 1/1 -> 0/0 as de novo (``--alt2ref-DNM``).
    ``include`` may use curly-bracket threshold ranges, for example ``GQ>{10,20,30}``. The report
    (``SN*``, ``FRM``, ``VAF`` and similar rows) is saved to ``output_path``.
    """
    options = [
        "-c", cli.text(csq_tag, "csq_tag"),
        "--max-len", str(cli.integer(max_length, "max_length", minimum=1)),
    ]
    if ped_path is not None:
        options += ["-p", str(cli.existing_file(ped_path, "ped_path"))]
    if h.flag(alt_to_ref_dnm, "alt_to_ref_dnm"):
        options.append("--alt2ref-DNM")
    return h.run_text_plugin(
        "bcftools_indel_stats", "indel-stats", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        header="none", bcftools_path=bcftools_path,
    )


def bcftools_vcf_to_table(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    hide: Sequence[str] | None = None,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Print variants as text tables (``bcftools +vcf2table``): variant, INFO, genotype and consequence tables.

    ``hide`` lists the tables to leave out: ``HOM_REF``/``RR``, ``HET``/``AR``, ``NO_CALL``/``MISSING``,
    ``CSQ``/``VEP``, ``SPLICEAI``, ``ANN``/``SNPEFF``, ``LOF``, ``VC``, ``INFO``, ``GT``, ``GTTYPES``,
    ``URL``. The text is saved to ``output_path``; ``variants`` counts the ``<<<`` headers.
    """
    options: list[str] = []
    if hide is not None:
        options += ["--hide", ",".join(h.text_list(list(hide), "hide"))]
    run = h.run_text_plugin(
        "bcftools_vcf_to_table", "vcf2table", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude),
        header="none", bcftools_path=bcftools_path,
    )
    text = h.read_text(Path(run["outputs"]["report"]))
    run["report"] = {"variants": sum(1 for line in text.splitlines() if line.startswith("<<<"))}
    run["text_preview"] = text[:2000]
    return run


def bcftools_variant_read_frequency_score(
    output_path: str | Path,
    *,
    alignment_list: str | Path | None = None,
    reference_fasta: str | Path | None = None,
    sites_file: str | Path | None = None,
    min_depth: int = 10,
    vaf_bins: int = 20,
    batch: str | None = None,
    merge_batch_files: Sequence[str | Path] | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variant read frequency score (site noisiness) from many BAM/CRAM files (``bcftools +vrfs``).

    Scoring mode: ``alignment_list`` is a text file with one BAM/CRAM path per line (``-a``),
    ``reference_fasta`` the reference (``-f``) and ``sites_file`` a tab-separated file of
    ``chrom, pos, ref, alt`` sites (``-s``); ``min_depth`` (``-d``) and ``vaf_bins`` (``-n``) tune the
    model. ``batch`` is ``I/N`` to run only the I-th of N batches (``-b``); the scores are written to
    ``output_path``. Merge mode: ``merge_batch_files`` lists batch outputs to merge (``-M``), and the
    scoring arguments must be left out.
    """
    target = cli.new_file(output_path, "output_path")
    if merge_batch_files is not None:
        if any(item is not None for item in (alignment_list, reference_fasta, sites_file, batch)):
            raise ConfigurationError("merge_batch_files cannot be combined with scoring arguments.")
        files = h.existing_files(list(merge_batch_files), "merge_batch_files")
        arguments = ["+vrfs", "-M", *(str(path) for path in files)]
    else:
        if alignment_list is None or reference_fasta is None or sites_file is None:
            raise ConfigurationError("alignment_list, reference_fasta and sites_file are required.")
        arguments = [
            "+vrfs",
            "-a", str(cli.existing_file(alignment_list, "alignment_list")),
            "-f", str(cli.existing_file(reference_fasta, "reference_fasta")),
            "-s", str(cli.existing_file(sites_file, "sites_file")),
            "-d", str(cli.integer(min_depth, "min_depth", minimum=0)),
            "-n", str(cli.integer(vaf_bins, "vaf_bins", minimum=2)),
        ]
        if batch is not None:
            if not isinstance(batch, str) or not re.match(r"^\d+/\d+$", batch):
                raise ConfigurationError("batch must look like 'I/N', for example '1/2'.")
            arguments += ["-b", batch]
    done = h.invoke("bcftools", bcftools_path, arguments, plugin=True, stdout_path=target)
    return h.result(
        "bcftools_variant_read_frequency_score", done, outputs={"report": target},
        report=h.tab_table(h.read_text(target), header="none"),
    )


# ---------------------------------------------------------------------------
# VCF-writing plugins
# ---------------------------------------------------------------------------


def bcftools_add_variantkey(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add the INFO fields ``VKX`` (hexadecimal VariantKey) and ``RSX`` (hexadecimal rsID number)
    (``bcftools +add-variantkey``).
    """
    return h.run_vcf_plugin(
        "bcftools_add_variantkey", "add-variantkey", vcf_path, output_path,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_variantkey_tables(
    vcf_path: str | Path,
    output_dir: str | Path,
    *,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate unsorted VariantKey lookup tables (``bcftools +variantkey-hex``).

    Writes ``vkrs.unsorted.hex`` (VariantKey to rsID), ``rsvk.unsorted.hex`` (rsID to VariantKey) and
    ``nrvk.unsorted.tsv`` (non-reversible VariantKeys) into the new folder ``output_dir``. Returns the
    files and the counters bcftools prints.
    """
    source = cli.existing_file(vcf_path, "vcf_path")
    folder = cli.new_directory(output_dir, "output_dir")
    try:
        done = h.invoke(
            "bcftools", bcftools_path, h.plugin_command("variantkey-hex", [source], separator=False),
            plugin=True, cwd=folder,
        )
    except BaseException:
        for leftover in folder.iterdir():
            leftover.unlink()
        folder.rmdir()
        raise
    counters = {
        key.strip(): int(value)
        for key, value in re.findall(r"^(.+?):\s+(\d+)\s*$", done.stdout + "\n" + done.stderr, re.MULTILINE)
    }
    return h.result(
        "bcftools_variantkey_tables", done, outputs={"output_dir": folder},
        files=h.directory_files(folder), counters=counters,
    )


def bcftools_contrast_groups(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    control_samples: Sequence[str],
    case_samples: Sequence[str],
    annotations: Sequence[str] | None = None,
    max_allele_freq: float | None = None,
    force_samples: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Association test and novel allele/genotype detection between two sample groups
    (``bcftools +contrast``).

    ``control_samples`` is the background group (``-0``), ``case_samples`` the group where novel alleles
    or genotypes are expected (``-1``). ``annotations`` chooses the INFO tags written, from ``PASSOC``
    (Fisher test, REF vs non-REF), ``FASSOC`` (non-REF fraction in controls and cases), ``NASSOC``
    (allele counts), ``NOVELAL`` and ``NOVELGT`` (default PASSOC, FASSOC, NOVELAL). ``max_allele_freq``
    (``-f``) tests enrichment of rare alleles (a value in (0, 1) is an ALT frequency, an integer an ALT
    count). ``force_samples`` continues when listed samples are missing from the VCF.
    """
    options = [
        "-0", ",".join(h.text_list(list(control_samples), "control_samples")),
        "-1", ",".join(h.text_list(list(case_samples), "case_samples")),
    ]
    if annotations is not None:
        items = h.text_list(list(annotations), "annotations")
        for item in items:
            cli.choice(item, _ANNOTATIONS, "annotations")
        options += ["-a", ",".join(items)]
    if max_allele_freq is not None:
        options += ["-f", _number_text(max_allele_freq, "max_allele_freq", minimum=0.0)]
    if h.flag(force_samples, "force_samples"):
        options.append("--force-samples")
    return h.run_vcf_plugin(
        "bcftools_contrast_groups", "contrast", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_fill_from_fasta(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    fasta_path: str | Path,
    column: str,
    header_lines_path: str | Path | None = None,
    include: str | None = None,
    exclude: str | None = None,
    replace_non_acgtn: bool = False,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fill the REF column or an INFO tag from a FASTA file (``bcftools +fill-from-fasta``).

    ``column`` is ``REF`` (repair REF alleles that are N or wrong) or an INFO tag such as ``AA`` for an
    ancestral allele; the FASTA must be indexed (``samtools faidx``). ``header_lines_path`` holds the
    ``##INFO`` line for a new tag (``-h``). ``include`` / ``exclude`` limit the annotated records.
    ``replace_non_acgtn`` replaces other characters with N (``-N``).
    """
    options = [
        "-c", cli.text(column, "column"),
        "-f", str(cli.existing_file(fasta_path, "fasta_path")),
    ]
    if header_lines_path is not None:
        options += ["-h", str(cli.existing_file(header_lines_path, "header_lines_path"))]
    if include is not None and exclude is not None:
        raise ConfigurationError("include and exclude are mutually exclusive.")
    if include is not None:
        options += ["-i", cli.text(include, "include")]
    if exclude is not None:
        options += ["-e", cli.text(exclude, "exclude")]
    if h.flag(replace_non_acgtn, "replace_non_acgtn"):
        options.append("-N")
    return h.run_vcf_plugin(
        "bcftools_fill_from_fasta", "fill-from-fasta", vcf_path, output_path, options=options,
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_fill_tags(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    tags: Sequence[str],
    populations_file: str | Path | None = None,
    drop_missing: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute INFO/FORMAT tags from genotypes and write the VCF (``bcftools +fill-tags``).

    ``tags`` selects what is filled (``-t``): ``AF``, ``AC``, ``AC_Hom``, ``AC_Het``, ``AC_Hemi``,
    ``AN``, ``ExcHet``, ``HWE``, ``MAF``, ``NS``, ``F_MISSING``, ``END``, ``TYPE``, ``VAF``, ``VAF1``,
    ``all``, or custom expressions such as ``DP:1=int(sum(FORMAT/DP))``. ``populations_file`` (``-S``)
    lists samples with comma-separated population names to compute the tags per population.
    ``drop_missing`` does not count half-missing genotypes as hemizygous (``-d``).
    """
    items = [cli.text(item, "tags") for item in (tags if not isinstance(tags, str) else [tags])]
    if not items:
        raise ConfigurationError("tags must not be empty.")
    options = ["-t", ",".join(items)]
    if populations_file is not None:
        options += ["-S", str(cli.existing_file(populations_file, "populations_file"))]
    if h.flag(drop_missing, "drop_missing"):
        options.append("-d")
    return h.run_vcf_plugin(
        "bcftools_fill_tags", "fill-tags", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude),
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_exclusive_genotype_sites(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    samples: Sequence[str],
    regions: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Output only the sites where the listed samples all share one genotype and no other sample has it
    (``bcftools +GTsubset``).
    """
    return h.run_vcf_plugin(
        "bcftools_exclusive_genotype_sites", "GTsubset", vcf_path, output_path,
        options=["-s", ",".join(h.text_list(list(samples), "samples"))],
        general=h.selection(regions=regions), output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_isec_genotypes(
    vcf_path: str | Path,
    other_vcf_path: str | Path,
    output_path: str | Path,
    *,
    regions: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set the genotypes of ``vcf_path`` that differ from ``other_vcf_path`` to missing (``bcftools +isecGT``).

    Both files must be bgzip-compressed or BCF and indexed; samples and sites are matched by name and
    position.
    """
    other = cli.existing_file(other_vcf_path, "other_vcf_path")
    return h.run_vcf_plugin(
        "bcftools_isec_genotypes", "isecGT", vcf_path, output_path,
        general=h.selection(regions=regions), separator=False, extra_inputs=[other],
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_prune_sites(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    statistic: str = "r2",
    max_value: float | None = None,
    annotate: Sequence[str] | None = None,
    window: str = "100kb",
    sites_per_window: int | None = None,
    sites_per_window_mode: str = "maxAF",
    af_tag: str | None = None,
    set_filter: str | None = None,
    keep_sites: bool = False,
    random_seed: int | None = None,
    randomize_missing: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prune or annotate sites by LD or site density in windows (``bcftools +prune``).

    ``statistic`` is ``r2`` (correlation), ``LD`` (Lewontin's D), ``RD`` (Ragsdale's D) or ``count``
    (number of variants in the window); ``max_value`` removes sites above it (``-m``), or, with
    ``set_filter`` (``-f``), marks them with that FILTER instead. ``annotate`` adds the biggest
    ``r2``/``LD``/``RD`` value and its position, or the ``count`` (``-a``). ``window`` is ``100kb``, ``1000``
    (sites) and similar (``-w``). ``sites_per_window`` keeps at most N sites per window (``-n``) chosen by
    ``sites_per_window_mode`` ``maxAF``, ``1st`` or ``rand``. ``af_tag`` (``--AF-tag``), ``keep_sites``
    (``-k``), ``random_seed`` and ``randomize_missing`` are the matching options. At least one of
    ``max_value``, ``sites_per_window`` or ``annotate`` is required.
    """
    kind = cli.choice(statistic, _PRUNE_STATS, "statistic")
    options: list[str] = ["-w", cli.text(window, "window")]
    if max_value is not None:
        options += ["-m", f"{kind}={_number_text(max_value, 'max_value', minimum=0.0)}"]
    if annotate is not None:
        items = h.text_list(list(annotate), "annotate")
        for item in items:
            cli.choice(item, _PRUNE_STATS, "annotate")
        options += ["-a", ",".join(items)]
    if sites_per_window is not None:
        options += ["-n", str(cli.integer(sites_per_window, "sites_per_window", minimum=1))]
        options += ["-N", cli.choice(sites_per_window_mode, _NSITES_MODES, "sites_per_window_mode")]
    if max_value is None and sites_per_window is None and annotate is None:
        raise ConfigurationError("Give max_value, sites_per_window or annotate.")
    if af_tag is not None:
        options += ["--AF-tag", cli.text(af_tag, "af_tag")]
    if set_filter is not None:
        options += ["-f", cli.text(set_filter, "set_filter")] if set_filter != "." else ["-f", "."]
    if h.flag(keep_sites, "keep_sites"):
        options.append("-k")
    if random_seed is not None:
        options += ["--random-seed", str(cli.integer(random_seed, "random_seed", minimum=0))]
    if h.flag(randomize_missing, "randomize_missing"):
        options.append("--randomize-missing")
    return h.run_vcf_plugin(
        "bcftools_prune_sites", "prune", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude), separator=False,
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_remove_overlaps(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    mark: str = "overlap",
    mark_tag: str | None = None,
    missing_value: str = "0",
    reverse: bool = False,
    output_type: str | None = None,
    regions: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove, list or mark overlapping and duplicate variants (``bcftools +remove-overlaps``).

    ``mark`` chooses the sites (``-m``): ``overlap`` (overlapping sites, default), ``dup`` (all duplicate
    sites) or ``min(QUAL)`` (the lowest-QUAL sites until the overlaps are resolved; ``missing_value``
    ``0`` or ``DP`` sets the QUAL of missing values). By default the chosen sites are removed; with
    ``mark_tag`` they are kept and flagged ``INFO/TAG`` instead (``-M``). ``reverse`` applies the
    reverse logic (keep duplicates instead of removing them). ``output_type`` ``t`` / ``tz`` writes a
    plain (compressed) list of the sites (chromosome, position) instead of a VCF.
    """
    if mark not in {"overlap", "dup", "min(QUAL)"}:
        raise ConfigurationError("mark must be overlap, dup or min(QUAL).")
    options = ["-m", mark, "--missing", cli.choice(missing_value, ("0", "DP"), "missing_value")]
    if mark_tag is not None:
        options += ["-M", cli.text(mark_tag, "mark_tag")]
    if h.flag(reverse, "reverse"):
        options.append("--reverse")
    if output_type in {"t", "tz"}:
        source = cli.existing_file(vcf_path, "vcf_path")
        target = cli.new_file(output_path, "output_path")
        arguments = h.plugin_command(
            "remove-overlaps", [source],
            general=[*h.selection(regions=regions), "-O", output_type, "-o", str(target)],
            options=options, separator=False,
        )
        done = h.invoke("bcftools", bcftools_path, arguments, plugin=True, created=[target])
        return h.result(
            "bcftools_remove_overlaps", done, outputs={"output": target}, output_type=output_type,
            messages=[line for line in done.stderr_tail.splitlines() if line.strip()][-4:],
        )
    return h.run_vcf_plugin(
        "bcftools_remove_overlaps", "remove-overlaps", vcf_path, output_path, options=options,
        general=h.selection(regions=regions), separator=False, output_type=output_type,
        bcftools_path=bcftools_path,
    )


def bcftools_scatter_vcf(
    vcf_path: str | Path,
    output_dir: str | Path,
    *,
    sites_per_chunk: int | None = None,
    regions_to_scatter: Sequence[str] | None = None,
    scatter_file: str | Path | None = None,
    extra_name: str | None = None,
    prefix: str | None = None,
    output_type: str = "v",
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a VCF/BCF into chunks of N sites or into region windows (``bcftools +scatter``).

    Exactly one of ``sites_per_chunk`` (``-n``), ``regions_to_scatter`` (``-s``, region strings, one
    output per region) and ``scatter_file`` (``-S``, regions with an optional second column naming the
    output) is required. ``extra_name`` writes the records outside the listed regions to a separate
    file of that name (``-x``); ``prefix`` is prepended to the output names (``-p``). The files are
    written into the new folder ``output_dir``; ``output_type`` is ``v``, ``z``, ``b`` or ``u``.
    """
    chosen = [value is not None for value in (sites_per_chunk, regions_to_scatter, scatter_file)]
    if sum(chosen) != 1:
        raise ConfigurationError("Give exactly one of sites_per_chunk, regions_to_scatter, scatter_file.")
    source = cli.existing_file(vcf_path, "vcf_path")
    options: list[str] = []
    if sites_per_chunk is not None:
        options += ["-n", str(cli.integer(sites_per_chunk, "sites_per_chunk", minimum=1))]
    if regions_to_scatter is not None:
        options += ["-s", ",".join(h.text_list(list(regions_to_scatter), "regions_to_scatter"))]
    if scatter_file is not None:
        options += ["-S", str(cli.existing_file(scatter_file, "scatter_file"))]
    if extra_name is not None:
        options += ["-x", cli.text(extra_name, "extra_name")]
    if prefix is not None:
        options += ["-p", cli.text(prefix, "prefix")]
    kind = cli.choice(output_type, h.OUTPUT_TYPES, "output_type")
    folder = cli.new_directory(output_dir, "output_dir")
    try:
        arguments = h.plugin_command(
            "scatter", [source], general=["-O", kind, "-o", str(folder)], options=options, separator=False
        )
        done = h.invoke("bcftools", bcftools_path, arguments, plugin=True, created=[])
    except BaseException:
        for leftover in folder.iterdir():
            leftover.unlink()
        folder.rmdir()
        raise
    return h.result(
        "bcftools_scatter_vcf", done, outputs={"output_dir": folder}, output_type=kind,
        files=h.directory_files(folder),
    )


def bcftools_set_genotypes(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    target_genotypes: str = ".",
    new_genotype: str,
    genotype_include: str | None = None,
    genotype_exclude: str | None = None,
    random_seed: int = 0,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set genotypes by rule (``bcftools +setGT``; also the ``+missing2ref`` plugin).

    ``target_genotypes`` (``-t``) picks the genotypes to change: ``./.`` (completely missing),
    ``./x`` (partly missing), ``.`` (partly or completely missing), ``a`` (all), ``q`` (the genotypes
    for which the expression ``genotype_include`` / ``genotype_exclude`` is true, plugin ``-i`` / ``-e``,
    for example ``GT="." && FMT/DP>0``; ``include`` / ``exclude`` select whole records instead), ``b`` (heterozygous genotypes failing a
    binomial test, for example ``b:AD<1e-3``) or ``r:FLOAT`` (a random proportion). ``new_genotype``
    (``-n``) is ``.`` (missing), ``0`` (reference), ``c:GT`` (custom genotype such as ``c:0/1``), ``i``
    (invert phase), ``m`` / ``M`` (minor / major allele), ``X`` (allele with bigger depth), ``p`` (phase)
    or ``u`` (unphase), and can combine ``0``/``m``/``M`` with ``p`` (``0p``). The ``+missing2ref``
    plugin is ``target_genotypes="."`` with ``new_genotype`` ``0`` (``0p`` for its ``-p``) or ``M`` (``Mp``
    for ``-m -p``).
    """
    if not isinstance(target_genotypes, str) or not _SET_TARGET.match(target_genotypes):
        raise ConfigurationError(
            "target_genotypes must be one of ./., ./x, ., a, q, b[:AD<p], r:FLOAT."
        )
    if not isinstance(new_genotype, str) or not _SET_NEW.match(new_genotype):
        raise ConfigurationError("new_genotype must be ., 0, c:GT, i, m, M, X, p or u (with an optional p).")
    options = [
        "-t", target_genotypes, "-n", new_genotype,
        "-s", str(cli.integer(random_seed, "random_seed", minimum=0)),
    ]
    if genotype_include is not None and genotype_exclude is not None:
        raise ConfigurationError("genotype_include and genotype_exclude are mutually exclusive.")
    if genotype_include is not None:
        options += ["-i", cli.text(genotype_include, "genotype_include")]
    if genotype_exclude is not None:
        options += ["-e", cli.text(genotype_exclude, "genotype_exclude")]
    if target_genotypes == "q" and genotype_include is None and genotype_exclude is None:
        raise ConfigurationError("target_genotypes 'q' needs genotype_include or genotype_exclude.")
    return h.run_vcf_plugin(
        "bcftools_set_genotypes", "setGT", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude),
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_split_vcf_by_sample(
    vcf_path: str | Path,
    output_dir: str | Path,
    *,
    samples_file: str | Path | None = None,
    groups_file: str | Path | None = None,
    keep_tags: Sequence[str] | None = None,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str = "v",
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a VCF/BCF by sample or by sample group (``bcftools +split``).

    By default one single-sample file is written per sample, named after the sample. ``samples_file``
    (``-S``) lists samples per output file (one line per file; optional second column renames, optional
    third column names the file); ``groups_file`` (``-G``) assigns each sample (second column renames, third
    column lists the files) to one or more outputs. ``keep_tags`` limits the tags kept, for example
    ``["INFO", "FMT/GT"]`` (``-k``). The files are written into the new folder ``output_dir``.
    """
    if samples_file is not None and groups_file is not None:
        raise ConfigurationError("samples_file and groups_file are mutually exclusive.")
    source = cli.existing_file(vcf_path, "vcf_path")
    options: list[str] = []
    if samples_file is not None:
        options += ["-S", str(cli.existing_file(samples_file, "samples_file"))]
    if groups_file is not None:
        options += ["-G", str(cli.existing_file(groups_file, "groups_file"))]
    if keep_tags is not None:
        options += ["-k", ",".join(h.text_list(list(keep_tags), "keep_tags"))]
    kind = cli.choice(output_type, h.OUTPUT_TYPES, "output_type")
    folder = cli.new_directory(output_dir, "output_dir")
    try:
        arguments = h.plugin_command(
            "split", [source],
            general=[*h.selection(regions=regions, include=include, exclude=exclude), "-O", kind, "-o", str(folder)],
            options=options, separator=False,
        )
        done = h.invoke("bcftools", bcftools_path, arguments, plugin=True)
    except BaseException:
        for leftover in folder.iterdir():
            leftover.unlink()
        folder.rmdir()
        raise
    return h.result(
        "bcftools_split_vcf_by_sample", done, outputs={"output_dir": folder}, output_type=kind,
        files=h.directory_files(folder),
    )


def bcftools_split_vep(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    annotation: str | None = None,
    columns: str | None = None,
    select: str | None = None,
    text_format: str | None = None,
    list_fields: bool = False,
    annotation_prefix: str | None = None,
    all_fields_delimiter: str | None = None,
    severity_file: str | Path | None = None,
    gene_list: str | Path | None = None,
    duplicate: bool = False,
    drop_sites: bool | None = None,
    allow_undefined_tags: bool = False,
    print_header: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Query or extract structured consequence annotations: VEP ``CSQ``, bcftools ``BCSQ``, SnpEff ``ANN``
    (``bcftools +split-vep``).

    ``annotation`` names the INFO tag (auto-detected: CSQ, BCSQ, ANN). ``list_fields`` (``-l``) lists the
    annotation fields instead of processing records. ``columns`` (``-c``) extracts fields by name or
    0-based index into new INFO tags (``Consequence,IMPACT,SYMBOL``, ``1-3``, ``-`` for all; add
    ``:Integer`` / ``:Float`` for types); ``select`` (``-s``) is ``TR:CSQ[:PRN]`` such as ``worst``,
    ``primary:missense+`` or ``all:any:worst``. With ``text_format`` (a ``bcftools query`` format such as
    ``%CHROM\\t%POS\\t%Consequence\\n``, ``-f``) the output is a text table (written to ``output_path``)
    instead of a VCF. ``annotation_prefix`` (``-p``), ``all_fields_delimiter`` (``-A``), ``severity_file``
    (``-S``), ``gene_list`` (``-g``), ``duplicate`` (``-d``, one line per transcript), ``drop_sites``
    (True ``-x`` / False ``-X``), ``allow_undefined_tags`` (``-u``) and ``print_header`` (``-H``) are the
    matching options.
    """
    options: list[str] = []
    if annotation is not None:
        options += ["-a", cli.text(annotation, "annotation")]
    if h.flag(list_fields, "list_fields"):
        options.append("-l")
    if columns is not None:
        options += ["-c", columns if columns == "-" else cli.text(columns, "columns")]
    if select is not None:
        options += ["-s", cli.text(select, "select")]
    if annotation_prefix is not None:
        options += ["-p", cli.text(annotation_prefix, "annotation_prefix")]
    if all_fields_delimiter is not None:
        options += ["-A", cli.text(all_fields_delimiter, "all_fields_delimiter")]
    if severity_file is not None:
        options += ["-S", str(cli.existing_file(severity_file, "severity_file"))]
    if gene_list is not None:
        options += ["-g", str(cli.existing_file(gene_list, "gene_list"))]
    for enabled, option, name in (
        (duplicate, "-d", "duplicate"),
        (allow_undefined_tags, "-u", "allow_undefined_tags"),
        (print_header, "-H", "print_header"),
    ):
        if h.flag(enabled, name):
            options.append(option)
    if drop_sites is not None:
        options.append("-x" if h.flag(drop_sites, "drop_sites") else "-X")
    general = h.selection(regions=regions, include=include, exclude=exclude)
    if text_format is not None or list_fields:
        if text_format is not None:
            options += ["-f", cli.text(text_format.replace("\\t", "\t").replace("\\n", "\n"), "text_format")]
        run = h.run_text_plugin(
            "bcftools_split_vep", "split-vep", vcf_path, output_path, options=options,
            general=general, separator=False, header="none", bcftools_path=bcftools_path,
        )
        if list_fields:
            fields = [
                {"index": row[0], "name": row[1]}
                for row in run["report"]["rows"]
                if len(row) >= 2 and isinstance(row[0], int)
            ]
            run["fields"] = fields
        return run
    return h.run_vcf_plugin(
        "bcftools_split_vep", "split-vep", vcf_path, output_path, options=options, general=general,
        separator=False, output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_convert_format_tags(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    conversion: str,
    replace_source: bool = False,
    gp_threshold: float = 0.1,
    skip_nalt: int = 0,
    defaults: str | None = None,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert between equivalent FORMAT tags (``bcftools +tag2tag``).

    ``conversion`` is one of ``GL-to-PL``, ``GL-to-GP``, ``GL-to-GT``, ``PL-to-GL``, ``PL-to-GP``,
    ``PL-to-GT``, ``GP-to-GL``, ``GP-to-PL``, ``GP-to-GT`` (the ``--ORI-to-NEW`` family), ``QR-QA-to-QS``,
    and the localized-allele conversions ``LXX-to-XX`` (``LAA,LPL,LAD`` to ``PL,AD``) or a single tag
    (``LPL-to-PL``, ``LAD-to-AD``). The opposite direction (``--XX-to-LXX``) is listed in the bcftools 1.24
    help but not implemented upstream ("todo"), so it is not offered. ``replace_source`` drops the source
    tag (``-r``). ``gp_threshold`` is the maximum error for a GP to GT hard call (``-t``). ``skip_nalt``
    leaves sites with more than that many ALT alleles unchanged in the localized conversions (``-s``,
    0 = all sites). ``defaults`` is
    ``AD:.,PL:.`` style values for missing localized fields (``-d``).
    """
    options = [f"--{cli.choice(conversion, _TAG2TAG, 'conversion')}"]
    if h.flag(replace_source, "replace_source"):
        options.append("-r")
    options += [
        "-t", _number_text(gp_threshold, "gp_threshold", minimum=0.0, maximum=1.0),
        "-s", str(cli.integer(skip_nalt, "skip_nalt", minimum=0)),
    ]
    if defaults is not None:
        options += ["-d", cli.text(defaults, "defaults")]
    return h.run_vcf_plugin(
        "bcftools_convert_format_tags", "tag2tag", vcf_path, output_path, options=options,
        general=h.selection(regions=regions, include=include, exclude=exclude),
        output_type=output_type, bcftools_path=bcftools_path,
    )


def bcftools_variant_distance(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    direction: str = "nearest",
    tag_name: str = "DIST",
    regions: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate sites with the distance to the nearest variant (``bcftools +variant-distance``).

    ``direction`` is ``nearest`` (default), ``fwd``, ``rev`` or ``both`` (two values: distance to the
    previous and to the next variant); ``tag_name`` is the INFO tag written (``-n``, default ``DIST``).
    """
    options = [
        "-d", cli.choice(direction, _DIRECTIONS, "direction"),
        "-n", cli.text(tag_name, "tag_name"),
    ]
    return h.run_vcf_plugin(
        "bcftools_variant_distance", "variant-distance", vcf_path, output_path, options=options,
        general=h.selection(regions=regions), separator=False, output_type=output_type,
        bcftools_path=bcftools_path,
    )


def bcftools_fixref(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    reference_fasta: str | Path,
    mode: str = "stats",
    dbsnp_vcf: str | Path | None = None,
    discard_unresolved: bool = False,
    tag_name: str | None = None,
    regions: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report or fix REF/ALT strand orientation and swaps against a reference (``bcftools +fixref``).

    ``mode`` (``-m``): ``stats`` (default) writes the strand and mismatch statistics report to
    ``output_path`` (the SC, ST, NS rows are parsed into ``report``); the conversion modes write a
    VCF/BCF: ``flip`` (flip or swap non-ambiguous SNPs), ``flip-all`` (also A/T and C/G sites), ``id``
    (use the ID column against ``dbsnp_vcf`` to decide), ``ref-alt`` (swap REF/ALT to match the
    reference), ``swap`` (swap columns, keep genotypes) and ``top`` (Illumina TOP strand to forward).
    ``discard_unresolved`` drops sites that cannot be repaired (``-d``); ``tag_name`` renames the
    ``INFO/FIXREF`` annotation (``-t``). Do not apply a conversion without understanding the strand
    convention of the data.
    """
    kind = cli.choice(mode, _FIXREF_MODES, "mode")
    options = ["-f", str(cli.existing_file(reference_fasta, "reference_fasta")), "-m", kind]
    if dbsnp_vcf is not None:
        options += ["-i", str(cli.existing_file(dbsnp_vcf, "dbsnp_vcf"))]
    elif kind == "id":
        raise ConfigurationError("mode 'id' needs dbsnp_vcf.")
    if h.flag(discard_unresolved, "discard_unresolved"):
        options.append("-d")
    if tag_name is not None:
        options += ["-t", cli.text(tag_name, "tag_name")]
    general = h.selection(regions=regions)
    if kind == "stats":
        source = cli.existing_file(vcf_path, "vcf_path")
        target = cli.new_file(output_path, "output_path")
        arguments = h.plugin_command("fixref", [source], general=general, options=options)
        done = h.invoke("bcftools", bcftools_path, arguments, plugin=True)
        with target.open("x", encoding="utf-8") as handle:  # the statistics are printed on standard error
            handle.write(done.stderr)
        table = h.tab_table(done.stderr, header="none")
        sections: dict[str, list[list[object]]] = {}
        for row in table["rows"]:
            sections.setdefault(str(row[0]), []).append(row[1:])
        return h.result(
            "bcftools_fixref", done, outputs={"report": target}, report=table, sections=sections, mode=kind
        )
    run = h.run_vcf_plugin(
        "bcftools_fixref", "fixref", vcf_path, output_path, options=options, general=general,
        output_type=output_type, bcftools_path=bcftools_path,
    )
    run["mode"] = kind
    return run


def bcftools_color_chromosomes(
    vcf_path: str | Path,
    output_dir: str | Path,
    *,
    trio: Sequence[str] | None = None,
    unrelated: Sequence[str] | None = None,
    prefix: str = "color-chrs",
    regions: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find shared chromosomal segments between a trio or two unrelated samples from phased genotypes
    (``bcftools +color-chrs``).

    ``trio`` is ``[mother, father, child]`` (``-t``) and ``unrelated`` two sample names (``-u``); give
    exactly one. The segments (``SG`` rows: chromosome, start, end and the haplotype labels) and the switch
    counts (``SW`` rows) are written to ``<prefix>.dat`` in the new folder ``output_dir`` (the format
    that bcftools' ``color-chrs.pl`` turns into a figure) and parsed into ``segments`` and ``switches``.
    Phased genotypes are required.
    """
    if (trio is None) == (unrelated is None):
        raise ConfigurationError("Give exactly one of trio and unrelated.")
    source = cli.existing_file(vcf_path, "vcf_path")
    stem = cli.text(prefix, "prefix")
    if "/" in stem:
        raise ConfigurationError("prefix must be a file name stem, not a path.")
    options = ["-p", stem]
    if trio is not None:
        names = h.text_list(list(trio), "trio")
        if len(names) != 3:
            raise ConfigurationError("trio must be [mother, father, child].")
        options += ["-t", ",".join(names)]
    if unrelated is not None:
        names = h.text_list(list(unrelated), "unrelated")
        if len(names) != 2:
            raise ConfigurationError("unrelated must name exactly two samples.")
        options += ["-u", ",".join(names)]
    folder = cli.new_directory(output_dir, "output_dir")
    try:
        arguments = h.plugin_command(
            "color-chrs", [source], general=h.selection(regions=regions), options=options
        )
        done = h.invoke("bcftools", bcftools_path, arguments, plugin=True, cwd=folder)
        data = folder / f"{stem}.dat"
        if not data.is_file():
            raise ConfigurationError("bcftools +color-chrs wrote no output (are the genotypes phased?).")
    except BaseException:
        for leftover in folder.iterdir():
            leftover.unlink()
        folder.rmdir()
        raise
    table = h.tab_table(h.read_text(data), header="none")
    return h.result(
        "bcftools_color_chromosomes", done, outputs={"data": data},
        segments=[row[1:] for row in table["rows"] if row and row[0] == "SG"],
        switches=[row[1:] for row in table["rows"] if row and row[0] == "SW"],
        files=h.directory_files(folder),
    )


__all__ = [
    "bcftools_ad_bias",
    "bcftools_add_variantkey",
    "bcftools_af_distribution",
    "bcftools_allele_length_counts",
    "bcftools_check_ploidy",
    "bcftools_check_sparsity",
    "bcftools_color_chromosomes",
    "bcftools_contrast_groups",
    "bcftools_convert_format_tags",
    "bcftools_exclusive_genotype_sites",
    "bcftools_fill_from_fasta",
    "bcftools_fill_tags",
    "bcftools_fixref",
    "bcftools_genotype_dosage",
    "bcftools_genotype_intersections",
    "bcftools_indel_stats",
    "bcftools_isec_genotypes",
    "bcftools_prune_sites",
    "bcftools_remove_overlaps",
    "bcftools_scatter_vcf",
    "bcftools_set_genotypes",
    "bcftools_split_vcf_by_sample",
    "bcftools_split_vep",
    "bcftools_variant_distance",
    "bcftools_variant_read_frequency_score",
    "bcftools_variantkey_tables",
    "bcftools_vcf_to_table",
]

_register_public_verbose(globals(), tuple(__all__))
