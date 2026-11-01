"""CRISPR screen analysis and off-target search facades (pack P16): MAGeCK, DrugZ, Cas-OFFinder.

``count_sgrnas_mageck``, ``rank_genes_mageck`` (RRA, optionally with copy-number correction or estimation),
``estimate_gene_essentiality_mageck`` (MLE), ``enrich_pathways_mageck``, ``plot_screen_results_mageck``,
``score_drug_gene_interactions_drugz`` (paired and unpaired) and ``search_offtargets_cas_offinder``. Programs run without a
shell through ``dnakit.bioinformatics._p16_cli`` and write to fresh output prefixes; results are JSON-serialisable.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.exceptions import ConfigurationError

_MAGECK = cli.Tool("mageck", "mageck", "Install MAGeCK (pip install mageck or conda install -c bioconda mageck) and put it on PATH or pass mageck_path.",
                   "0.5.9.5", "https://doi.org/10.1186/s13059-014-0554-4", ("--version",), None)
_DRUGZ = cli.Tool("drugz", sys.executable, "Download drugz.py from https://github.com/hart-lab/drugz (needs pandas, scipy and six) and pass drugz_script.",
                  "1.1.0.2", "https://doi.org/10.1186/s13073-019-0665-3", ("--version",), None)
_CAS = cli.Tool("cas-offinder", "cas-offinder", "Install Cas-OFFinder (conda install -c bioconda cas-offinder) with an OpenCL runtime (pocl) and put it on PATH or pass cas_offinder_path.",
                "2.4.1", "https://doi.org/10.1093/bioinformatics/btu048", ("--version",), r"v(\d+(?:\.\d+)+)")
_NORM = ("none", "median", "total", "control")


def _prefix(value: str, suffixes: Sequence[str]) -> Path:
    prefix = Path(value).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder of the output prefix does not exist: {prefix.parent}")
    for suffix in suffixes:
        if Path(str(prefix) + suffix).exists():
            raise ConfigurationError(f"Output already exists and will not be overwritten: {prefix}{suffix}")
    return prefix


def _labels(values: Sequence[str] | str, name: str) -> str:
    items = [values] if isinstance(values, str) else list(values)
    if not items:
        raise ConfigurationError(f"{name} must not be empty.")
    return ",".join(cli.text(item, name) for item in items)


def count_sgrnas_mageck(
    library: str,
    output_prefix: str,
    *,
    fastq: Sequence[str] | None = None,
    sample_labels: Sequence[str] | None = None,
    alignment_files: Sequence[str] | None = None,
    normalization: str = "median",
    trim_5: int | str | None = None,
    sgrna_length: int | None = None,
    reverse_complement: bool = False,
    count_n: bool = False,
    mageck_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count sgRNA reads from FASTQ (or aligned reads) against an sgRNA library (``mageck count``).

    ``library`` is the sgRNA table (id, sequence, gene). ``fastq`` is one FASTQ file per sample (gzip allowed), in the
    order of ``sample_labels``; ``alignment_files`` can be given instead of ``fastq`` for SAM/BAM input. ``normalization`` is none, median, total or control. Writes
    ``<prefix>.count.txt``, ``<prefix>.count_normalized.txt`` and ``<prefix>.countsummary.txt`` (parsed into ``summary``).
    """
    lib = cli.existing_file(library, "library")
    prefix = _prefix(output_prefix, [".count.txt", ".countsummary.txt"])
    inputs = fastq if fastq is not None else alignment_files
    if inputs is None or (fastq is not None and alignment_files is not None):
        raise ConfigurationError("Give exactly one of fastq and alignment_files.")
    files = cli.existing_files(list(inputs), "fastq" if fastq is not None else "alignment_files")
    arguments = ["count", "--list-seq", str(lib), "--output-prefix", str(prefix), "--fastq", *map(str, files),
                 "--norm-method", cli.choice(normalization, _NORM, "normalization")]
    if sample_labels is not None:
        arguments += ["--sample-label", _labels(sample_labels, "sample_labels")]
    if trim_5 is not None:
        arguments += ["--trim-5", str(trim_5) if isinstance(trim_5, str) else str(cli.integer(trim_5, "trim_5", minimum=0))]
    if sgrna_length is not None:
        arguments += ["--sgrna-len", str(cli.integer(sgrna_length, "sgrna_length", minimum=1))]
    if cli.flag(reverse_complement, "reverse_complement"):
        arguments.append("--reverse-complement")
    if cli.flag(count_n, "count_n"):
        arguments.append("--count-n")
    expected = [Path(str(prefix) + ".count.txt"), Path(str(prefix) + ".countsummary.txt")]
    record = cli.invoke(_MAGECK, arguments, cwd=prefix.parent, expected=expected, explicit=mageck_path)
    return cli.finish(record, "count_sgrnas_mageck", summary=cli.read_table(expected[1], text_columns=frozenset({"Label", "File"})),
                      count_table=str(expected[0]))


def rank_genes_mageck(
    count_table: str,
    output_prefix: str,
    *,
    treatment: Sequence[str] | str,
    control: Sequence[str] | str | None = None,
    day0_label: str | None = None,
    paired: bool = False,
    normalization: str = "median",
    adjust_method: str = "fdr",
    gene_lfc_method: str | None = None,
    cnv_normalization: str | None = None,
    cell_line: str | None = None,
    cnv_estimation: str | None = None,
    control_sgrna: str | None = None,
    remove_zero: str | None = None,
    write_normalized_counts: bool = False,
    mageck_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene ranking of a count table by robust rank aggregation (``mageck test``), negative and positive selection.

    ``treatment``/``control`` are sample labels (or indices) of the count table; ``day0_label`` compares to day 0 instead of
    ``control``. ``cnv_normalization`` is a copy-number matrix file (tab-separated, first column ``SYMBOL``, one column per cell line) and ``cell_line`` the column to
    use (``--cnv-norm``/``--cell-line``: copy-number-corrected sgRNA scores); ``cnv_estimation`` is a BED file of genome
    regions for ``--cnv-est`` (estimate copy-number profiles from the screen). Writes ``<prefix>.gene_summary.txt`` and
    ``<prefix>.sgrna_summary.txt``; the top 20 gene rows are in ``genes``.
    """
    table = cli.existing_file(count_table, "count_table")
    prefix = _prefix(output_prefix, [".gene_summary.txt", ".sgrna_summary.txt"])
    arguments = ["test", "--count-table", str(table), "--output-prefix", str(prefix), "--treatment-id", _labels(treatment, "treatment"),
                 "--norm-method", cli.choice(normalization, _NORM, "normalization"),
                 "--adjust-method", cli.choice(adjust_method, ("fdr", "holm", "pounds"), "adjust_method")]
    if control is not None:
        arguments += ["--control-id", _labels(control, "control")]
    if day0_label is not None:
        arguments += ["--day0-label", cli.text(day0_label, "day0_label")]
    if control is None and day0_label is None:
        raise ConfigurationError("Give control sample labels or day0_label.")
    if cli.flag(paired, "paired"):
        arguments.append("--paired")
    if gene_lfc_method is not None:
        arguments += ["--gene-lfc-method", cli.choice(gene_lfc_method, ("median", "alphamedian", "mean", "alphamean", "secondbest"), "gene_lfc_method")]
    if cnv_normalization is not None:
        arguments += ["--cnv-norm", str(cli.existing_file(cnv_normalization, "cnv_normalization"))]
        if cell_line is not None:
            arguments += ["--cell-line", cli.text(cell_line, "cell_line")]
    if cnv_estimation is not None:
        arguments += ["--cnv-est", str(cli.existing_file(cnv_estimation, "cnv_estimation"))]
    if control_sgrna is not None:
        arguments += ["--control-sgrna", str(cli.existing_file(control_sgrna, "control_sgrna"))]
    if remove_zero is not None:
        arguments += ["--remove-zero", cli.choice(remove_zero, ("none", "control", "treatment", "both", "any"), "remove_zero")]
    if cli.flag(write_normalized_counts, "write_normalized_counts"):
        arguments.append("--normcounts-to-file")
    expected = [Path(str(prefix) + ".gene_summary.txt"), Path(str(prefix) + ".sgrna_summary.txt")]
    record = cli.invoke(_MAGECK, arguments, cwd=prefix.parent, expected=expected, explicit=mageck_path)
    return cli.finish(record, "rank_genes_mageck", genes=cli.read_table(expected[0], limit=20, text_columns=frozenset({"id"})),
                      gene_summary=str(expected[0]), sgrna_summary=str(expected[1]))


def estimate_gene_essentiality_mageck(
    count_table: str,
    design_matrix: str,
    output_prefix: str,
    *,
    beta_labels: Sequence[str] | None = None,
    include_samples: Sequence[str] | None = None,
    normalization: str = "median",
    cnv_normalization: str | None = None,
    cell_line: str | None = None,
    cnv_estimation: str | None = None,
    control_sgrna: str | None = None,
    permutation_rounds: int | None = None,
    remove_outliers: bool = False,
    threads: int = 1,
    mageck_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Maximum-likelihood gene essentiality with a design matrix (``mageck mle``, EM with permutation p-values).

    ``design_matrix`` is the sample-by-condition table (first column sample label, ``baseline`` column of 1). The
    ``cnv_*`` parameters correct (``--cnv-norm`` with ``--cell-line``) or estimate (``--cnv-est``) copy-number effects as in
    :func:`rank_genes_mageck`. Writes ``<prefix>.gene_summary.txt`` and ``<prefix>.sgrna_summary.txt``.
    """
    table = cli.existing_file(count_table, "count_table")
    design = cli.existing_file(design_matrix, "design_matrix")
    prefix = _prefix(output_prefix, [".gene_summary.txt", ".sgrna_summary.txt"])
    arguments = ["mle", "--count-table", str(table), "--design-matrix", str(design), "--output-prefix", str(prefix),
                 "--norm-method", cli.choice(normalization, _NORM, "normalization"), "--threads", str(cli.threads_value(threads))]
    if beta_labels is not None:
        arguments += ["--beta-labels", _labels(beta_labels, "beta_labels")]
    if include_samples is not None:
        arguments += ["--include-samples", _labels(include_samples, "include_samples")]
    if cnv_normalization is not None:
        arguments += ["--cnv-norm", str(cli.existing_file(cnv_normalization, "cnv_normalization"))]
        if cell_line is not None:
            arguments += ["--cell-line", cli.text(cell_line, "cell_line")]
    if cnv_estimation is not None:
        arguments += ["--cnv-est", str(cli.existing_file(cnv_estimation, "cnv_estimation"))]
    if control_sgrna is not None:
        arguments += ["--control-sgrna", str(cli.existing_file(control_sgrna, "control_sgrna"))]
    if permutation_rounds is not None:
        arguments += ["--permutation-round", str(cli.integer(permutation_rounds, "permutation_rounds", minimum=1))]
    if cli.flag(remove_outliers, "remove_outliers"):
        arguments.append("--remove-outliers")
    expected = [Path(str(prefix) + ".gene_summary.txt"), Path(str(prefix) + ".sgrna_summary.txt")]
    record = cli.invoke(_MAGECK, arguments, cwd=prefix.parent, expected=expected, explicit=mageck_path)
    return cli.finish(record, "estimate_gene_essentiality_mageck",
                      genes=cli.read_table(expected[0], limit=20, text_columns=frozenset({"Gene"})), gene_summary=str(expected[0]))


def enrich_pathways_mageck(
    gene_ranking: str,
    gene_sets: str,
    output_prefix: str,
    *,
    method: str = "gsea",
    single_ranking: bool = False,
    ranking_column: str | int | None = None,
    sort_criteria: str | None = None,
    permutations: int | None = None,
    mageck_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pathway/gene-set enrichment of a ranked gene list (``mageck pathway``) with GSEA or RRA.

    ``gene_ranking`` is a ``mageck test``/``mle`` gene summary (or, with ``single_ranking``, a two-column gene ranking);
    ``gene_sets`` is a GMT file. Writes ``<prefix>.pathway_summary.txt`` (top rows in ``pathways``).
    """
    ranking = cli.existing_file(gene_ranking, "gene_ranking")
    gmt = cli.existing_file(gene_sets, "gene_sets")
    prefix = _prefix(output_prefix, [".pathway_summary.txt"])
    arguments = ["pathway", "--gene-ranking", str(ranking), "--gmt-file", str(gmt), "--output-prefix", str(prefix),
                 "--method", cli.choice(method, ("gsea", "rra"), "method")]
    if cli.flag(single_ranking, "single_ranking"):
        arguments.append("--single-ranking")
    if ranking_column is not None:
        arguments += ["--ranking-column", str(ranking_column) if isinstance(ranking_column, int) else cli.text(ranking_column, "ranking_column")]
    if sort_criteria is not None:
        arguments += ["--sort-criteria", cli.choice(sort_criteria, ("neg", "pos"), "sort_criteria")]
    if permutations is not None:
        arguments += ["--permutation", str(cli.integer(permutations, "permutations", minimum=1))]
    summary = Path(str(prefix) + ".pathway_summary.txt")
    record = cli.invoke(_MAGECK, arguments, cwd=prefix.parent, expected=[summary], explicit=mageck_path)
    return cli.finish(record, "enrich_pathways_mageck", pathways=cli.read_table(summary, limit=20, text_columns=frozenset({"Pathway"})))


def plot_screen_results_mageck(
    count_table: str,
    gene_summary: str,
    output_prefix: str,
    *,
    genes: Sequence[str],
    samples: Sequence[str] | None = None,
    normalization: str = "median",
    mageck_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-gene sgRNA count plots for the listed genes (``mageck plot``; MAGeCK runs R, so ``Rscript`` is required).

    Writes ``<prefix>_sgrnacounts.pdf`` and the data table ``<prefix>_summary.txt`` when the program produces them; the
    files actually written are listed in ``outputs``.
    """
    table = cli.existing_file(count_table, "count_table")
    summary = cli.existing_file(gene_summary, "gene_summary")
    prefix = _prefix(output_prefix, ["_sgrnacounts.pdf", ".R"])
    arguments = ["plot", "--count-table", str(table), "--gene-summary", str(summary), "--output-prefix", str(prefix),
                 "--genes", _labels(genes, "genes"), "--norm-method", cli.choice(normalization, _NORM, "normalization")]
    if samples is not None:
        arguments += ["--samples", _labels(samples, "samples")]
    record = cli.invoke(_MAGECK, arguments, cwd=prefix.parent, explicit=mageck_path)
    produced = sorted(str(p) for p in prefix.parent.glob(prefix.name + "*") if p.is_file())
    if not any(path.endswith(".pdf") for path in produced):
        raise cli.BackendExecutionError("mageck plot produced no PDF (is Rscript installed?).", context={"log_tail": record["log_tail"]})
    record["outputs"] = produced
    return cli.finish(record, "plot_screen_results_mageck", genes=list(genes))


def score_drug_gene_interactions_drugz(
    count_table: str,
    output_path: str,
    *,
    control_samples: Sequence[str],
    drug_samples: Sequence[str],
    drugz_script: str,
    unpaired: bool = False,
    fold_change_path: str | None = None,
    remove_genes: Sequence[str] | None = None,
    pseudocount: int = 5,
    min_observations: int = 1,
    half_window_size: int = 500,
    drugz_python: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score drug-gene interactions in a paired sgRNA count screen with DrugZ (``drugz.py``).

    ``count_table`` is the tab-separated sgRNA read count table (first column ``GENE_CLONE``, second the gene). With
    ``unpaired=False`` control and drug samples are compared pairwise in the order given; ``unpaired=True`` (``-unpaired``)
    compares the mean fold change of the drug against the control group. ``drugz_script`` is the downloaded ``drugz.py`` and
    ``drugz_python`` the interpreter that has pandas, scipy and six (default: the current one). The gene table is written to
    ``output_path`` (columns normZ, pval/rank/fdr for synergy and suppression), the guide fold changes to
    ``fold_change_path`` when given.
    """
    table = cli.existing_file(count_table, "count_table")
    script = cli.existing_file(drugz_script, "drugz_script")
    target = cli.new_file(output_path, "output_path")
    controls, drugs = _labels(control_samples, "control_samples"), _labels(drug_samples, "drug_samples")
    arguments = [str(script), "-i", str(table), "-o", str(target), "-c", controls, "-x", drugs,
                 "-p", str(cli.integer(pseudocount, "pseudocount", minimum=0)), "--minobs", str(cli.integer(min_observations, "min_observations", minimum=1)),
                 "--half_window_size", str(cli.integer(half_window_size, "half_window_size", minimum=1)), "-q"]
    expected = [target]
    if fold_change_path is not None:
        fold = cli.new_file(fold_change_path, "fold_change_path")
        arguments += ["-f", str(fold)]
        expected.append(fold)
    if remove_genes:
        arguments += ["-r", _labels(remove_genes, "remove_genes")]
    if cli.flag(unpaired, "unpaired"):
        arguments.append("-unpaired")
    record = cli.invoke(_DRUGZ, arguments, cwd=target.parent, expected=expected, explicit=drugz_python or sys.executable)
    rows = cli.read_table(target, limit=20, text_columns=frozenset({"GENE"}))
    record["version"] = "drugz.py (script run by " + Path(drugz_python or sys.executable).name + ")"
    return cli.finish(record, "score_drug_gene_interactions_drugz", unpaired=unpaired, genes=rows)


def search_offtargets_cas_offinder(
    genome: str,
    guides: Mapping[str, int],
    output_path: str,
    *,
    pattern: str = "NNNNNNNNNNNNNNNNNNNNNRG",
    device: str = "C",
    cas_offinder_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a genome for CRISPR off-target sites tolerant to mismatches and the PAM (``cas-offinder``).

    ``genome`` is a FASTA file or a folder of FASTA/2bit files. ``pattern`` is the IUPAC guide-plus-PAM pattern (default
    20 nt + NRG, i.e. SpCas9 NGG/NAG); ``guides`` maps each guide (same length as ``pattern``, PAM positions as N) to the
    maximum number of mismatches. ``device`` is ``C`` (CPU), ``G`` (GPU) or ``A`` (accelerators), optionally followed by
    device ids. Output rows: guide, chromosome, position (0-based), site, strand, mismatches (parsed into ``sites``,
    first 1,000). The Cas-OFFinder input file is written next to the output (``<output_path>.input.txt``).
    """
    source = cli.existing_path(genome, "genome")
    target = cli.new_file(output_path, "output_path")
    spec = cli.text(pattern, "pattern").upper()
    if not isinstance(guides, Mapping) or not guides:
        raise ConfigurationError("guides must be a non-empty mapping of guide sequence to maximum mismatches.")
    lines = [str(source), spec]
    for sequence, mismatches in guides.items():
        guide = cli.text(sequence, "guides").upper()
        if len(guide) != len(spec):
            raise ConfigurationError(f"Guide {guide!r} must have the same length as pattern ({len(spec)}).")
        lines.append(f"{guide} {cli.integer(mismatches, 'mismatches', minimum=0, maximum=len(spec))}")
    if not (isinstance(device, str) and device[:1] in {"C", "G", "A"}):
        raise ConfigurationError("device must start with C, G or A.")
    spec_file = target.with_name(target.name + ".input.txt")
    spec_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    record = cli.invoke(_CAS, [str(spec_file), device, str(target)], cwd=target.parent, expected=[target], explicit=cas_offinder_path)
    sites = cli.read_table(target, header=False, names=["guide", "chrom", "position", "site", "strand", "mismatches"],
                           text_columns=frozenset({"guide", "chrom", "site", "strand"}), limit=1000)
    return cli.finish(record, "search_offtargets_cas_offinder", site_count=cli.count_lines(target), sites=sites,
                      input_file=str(spec_file))


__all__ = [
    "count_sgrnas_mageck",
    "enrich_pathways_mageck",
    "estimate_gene_essentiality_mageck",
    "plot_screen_results_mageck",
    "score_drug_gene_interactions_drugz",
    "search_offtargets_cas_offinder",
    "rank_genes_mageck",
]

_register_public_verbose(globals(), tuple(__all__))
