"""QUAST command-line facades: reference-based assembly evaluation, gene and rRNA prediction, read mapping and MetaQUAST.

QUAST (5.x, tested against 5.2.0) evaluates genome assemblies. The user installs it (``conda install -c bioconda quast`` or
``pip install quast``) and puts ``quast.py`` (and ``metaquast.py``) on PATH or passes ``quast_path``. Each function runs the
program with a fixed argument whitelist and no shell, writes into a new output folder (created by QUAST) and returns a
JSON-serialisable dictionary: status, backend, version, command, output files, the parsed ``report`` (metric to value per
assembly) and a bounded log tail. Plots, HTML reports and Icarus viewers are not built unless ``build_reports=True``.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="quast",
    executable="quast.py",
    hint="Install QUAST 5.x (`conda install -c bioconda quast` or `pip install quast`) and put quast.py and metaquast.py on PATH or pass quast_path.",
    tested_version="5.2.0",
    citation="https://doi.org/10.1093/bioinformatics/btt086",
    version_arguments=("--version",),
    version_pattern=r"QUAST v?(\d+\.\d+\.\d+)",
)
_META_TOOL = _cli.Tool(
    name="metaquast",
    executable="metaquast.py",
    hint=_TOOL.hint,
    tested_version="5.2.0",
    citation="https://doi.org/10.1093/bioinformatics/btv697",
    version_arguments=("--version",),
    version_pattern=r"QUAST v?(\d+\.\d+\.\d+)",
)
_GENE_FINDERS = ("genemarks", "mgm", "glimmer")


def _output_directory(value: str, name: str = "output_directory") -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def _report(directory: Path) -> dict[str, dict[str, Any]]:
    """The ``report.tsv`` of a QUAST run as ``{metric: {assembly: value}}``."""
    path = directory / "report.tsv"
    if not path.is_file():
        return {}
    lines = [
        line.split("\t")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line
    ]
    assemblies = lines[0][1:]
    return {
        row[0]: {
            name: _cli.convert_number(value)
            for name, value in zip(assemblies, row[1:], strict=False)
        }
        for row in lines[1:]
    }


def _common(
    contigs: Sequence[str],
    threads: int,
    min_contig: int,
    eukaryote: bool,
    large: bool,
    labels: Sequence[str] | None,
    build_reports: bool,
) -> tuple[list[str], list[str]]:
    files = [str(path) for path in _cli.existing_files(contigs, "contigs")]
    arguments = [
        "-t",
        str(_cli.threads_value(threads)),
        "-m",
        str(_cli.integer(min_contig, "min_contig", minimum=0)),
    ]
    for enabled, option in ((eukaryote, "--eukaryote"), (large, "--large")):
        if _cli.flag(enabled, option):
            arguments.append(option)
    if labels is not None:
        if (
            isinstance(labels, str)
            or len(labels) != len(files)
            or any("," in str(item) for item in labels)
        ):
            raise ConfigurationError(
                "labels must be a list with one comma-free label per contig file."
            )
        arguments += ["--labels", ", ".join(_cli.text(item, "labels") for item in labels)]
    if not _cli.flag(build_reports, "build_reports"):
        arguments += ["--no-plots", "--no-html", "--no-icarus"]
    return arguments, files


def _run(
    tool: _cli.Tool,
    arguments: Sequence[str],
    files: Sequence[str],
    out: Path,
    quast_path: str | None,
    timeout_seconds: float,
    operation: str,
    **extra: Any,
) -> dict[str, Any]:
    record = _cli.invoke(
        tool,
        ["-o", str(out), *arguments, *files],
        cwd=out.parent,
        expected=[out / "report.tsv"],
        explicit=quast_path,
        timeout_seconds=timeout_seconds,
    )
    record["outputs"] = sorted(str(path) for path in out.glob("*") if path.is_file())
    return _cli.finish(record, operation, output_directory=str(out), report=_report(out), **extra)


def quast_assess_assembly(
    contigs: Sequence[str],
    output_directory: str,
    *,
    reference: str | None = None,
    features: str | None = None,
    operons: str | None = None,
    k_mer_stats: bool = False,
    k_mer_size: int = 101,
    upper_bound_assembly: bool = False,
    upper_bound_reads: Sequence[str] | None = None,
    sv_bedpe: str | None = None,
    circos: bool = False,
    min_alignment: int = 65,
    min_identity: float = 95.0,
    extensive_misassembly_size: int = 1000,
    eukaryote: bool = False,
    large: bool = False,
    min_contig: int = 500,
    labels: Sequence[str] | None = None,
    build_reports: bool = False,
    threads: int = 1,
    quast_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate assemblies against a reference: genome fraction, NGAx, misassemblies, mismatches, indels, genes and operons (``quast.py -r``).

    ``contigs`` are one or more assembly FASTA files compared in one run. Without ``reference`` QUAST reports only reference-free
    statistics (N50, L50, GC...). With ``reference`` the contigs are aligned (minimap2) and QUAST reports genome fraction,
    duplication ratio, NGA50/LGA50, misassemblies (``extensive_misassembly_size`` is its ``-x``), mismatches and indels per 100 kb
    (``min_alignment`` and ``min_identity`` are ``-i`` and ``--min-identity``). ``features`` (GFF, BED, NCBI or TXT; ``-g``,
    optionally ``type:file`` for one GFF feature type) adds the number of completely and partially covered genes, and ``operons``
    (``--operons``) the same for operons. ``k_mer_stats`` (``-k``, size ``k_mer_size``) adds k-mer-based completeness and
    misjoins. ``upper_bound_assembly`` simulates the best possible assembly from the reference and the reads in
    ``upper_bound_reads`` (``--pe1``/``--pe2`` pair or one long-read file, see ``quast_map_reads`` for the read options);
    ``sv_bedpe`` compares the detected structural variants with a BEDPE list; ``circos`` draws a Circos plot (needs
    ``build_reports``). Result keys: ``report`` (``{metric: {assembly: value}}``) and the files of the run.
    """
    out = _output_directory(output_directory)
    arguments, files = _common(
        contigs, threads, min_contig, eukaryote, large, labels, build_reports
    )
    if reference is not None:
        arguments += ["-r", str(_cli.existing_file(reference, "reference"))]
    elif any(item is not None for item in (features, operons, sv_bedpe)) or upper_bound_assembly:
        raise ConfigurationError(
            "features, operons, sv_bedpe and upper_bound_assembly need a reference."
        )
    if features is not None:
        arguments += ["-g", str(_cli.existing_file(features, "features"))]
    if operons is not None:
        arguments += ["--operons", str(_cli.existing_file(operons, "operons"))]
    if _cli.flag(k_mer_stats, "k_mer_stats"):
        arguments += ["-k", "--k-mer-size", str(_cli.integer(k_mer_size, "k_mer_size", minimum=3))]
    if _cli.flag(upper_bound_assembly, "upper_bound_assembly"):
        arguments.append("--upper-bound-assembly")
        if upper_bound_reads is None:
            raise ConfigurationError(
                "upper_bound_assembly needs upper_bound_reads (one or two FASTQ files)."
            )
        reads = [str(path) for path in _cli.existing_files(upper_bound_reads, "upper_bound_reads")]
        if len(reads) == 2:
            arguments += ["--pe1", reads[0], "--pe2", reads[1]]
        elif len(reads) == 1:
            arguments += ["--single", reads[0]]
        else:
            raise ConfigurationError(
                "upper_bound_reads must be one file (single reads) or two files (a pair)."
            )
    if sv_bedpe is not None:
        arguments += ["--sv-bedpe", str(_cli.existing_file(sv_bedpe, "sv_bedpe"))]
    if _cli.flag(circos, "circos"):
        arguments.append("--circos")
    arguments += [
        "-i",
        str(_cli.integer(min_alignment, "min_alignment", minimum=1)),
        "--min-identity",
        format(_cli.number(min_identity, "min_identity", minimum=80.0, maximum=100.0), ".12g"),
        "-x",
        str(_cli.integer(extensive_misassembly_size, "extensive_misassembly_size", minimum=1)),
    ]
    return _run(
        _TOOL,
        arguments,
        files,
        out,
        quast_path,
        timeout_seconds,
        "reference_based_assembly_evaluation",
    )


def quast_predict_genes(
    contigs: Sequence[str],
    output_directory: str,
    *,
    gene_finder: str | None = "genemarks",
    rna_finding: bool = False,
    gene_thresholds: Sequence[int] = (0, 300, 1500, 3000),
    reference: str | None = None,
    eukaryote: bool = False,
    min_contig: int = 500,
    threads: int = 1,
    quast_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict genes and rRNA genes on assembly contigs and report their counts (``quast.py -f`` and ``--rna-finding``).

    ``gene_finder`` is ``'genemarks'`` (GeneMarkS for prokaryotes, GeneMark-ES with ``eukaryote=True``; the default finder of
    ``-f``), ``'mgm'`` (MetaGeneMark, ``--mgm``) or ``'glimmer'`` (GlimmerHMM, ``--glimmer``); ``None`` skips gene finding.
    ``rna_finding`` predicts ribosomal RNA genes with Barrnap (``--rna-finding``). ``gene_thresholds`` are the gene length
    cut-offs of the report (``--gene-thresholds``). The report rows ``# predicted genes (unique)`` and ``# predicted
    rRNA genes`` are in ``report``; the predicted gene files are in ``predicted_genes/`` of the output folder.
    """
    out = _output_directory(output_directory)
    arguments, files = _common(contigs, threads, min_contig, eukaryote, False, None, False)
    if reference is not None:
        arguments += ["-r", str(_cli.existing_file(reference, "reference"))]
    if gene_finder is not None:
        chosen = _cli.choice(gene_finder, _GENE_FINDERS, "gene_finder")
        arguments.append("-f")
        if chosen == "mgm":
            arguments.append("--mgm")
        elif chosen == "glimmer":
            arguments.append("--glimmer")
        if isinstance(gene_thresholds, str) or not gene_thresholds:
            raise ConfigurationError("gene_thresholds must be a non-empty list of lengths.")
        arguments += [
            "--gene-thresholds",
            ",".join(
                str(_cli.integer(item, "gene_thresholds", minimum=0)) for item in gene_thresholds
            ),
        ]
    if _cli.flag(rna_finding, "rna_finding"):
        arguments.append("--rna-finding")
    if gene_finder is None and not rna_finding:
        raise ConfigurationError("give gene_finder or rna_finding.")
    return _run(
        _TOOL,
        arguments,
        files,
        out,
        quast_path,
        timeout_seconds,
        "predict_genes_on_contigs",
        gene_finder=gene_finder,
        rna_finding=rna_finding,
    )


def quast_map_reads(
    contigs: Sequence[str],
    output_directory: str,
    *,
    reference: str | None = None,
    pe1: str | None = None,
    pe2: str | None = None,
    single: str | None = None,
    pacbio: str | None = None,
    nanopore: str | None = None,
    reference_bam: str | None = None,
    contig_bam_files: Sequence[str] | None = None,
    detect_variants: bool = True,
    min_contig: int = 500,
    threads: int = 1,
    quast_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map reads to the assembly and the reference: coverage, read-based SNPs, indels and structural variants (QUAST read options).

    Give reads as a paired-end pair (``pe1`` and ``pe2``), ``single``, ``pacbio`` or ``nanopore`` FASTQ files, or give the
    alignments directly: ``reference_bam`` (reads aligned to the reference) and ``contig_bam_files`` (one BAM per contig file,
    comma-joined for ``--bam``). QUAST aligns the reads, reports mapped reads, coverage depth and, with a ``reference`` and
    ``detect_variants`` (``False`` is ``--no-snps --no-sv``), SNPs, indels and structural variants. The numbers are in
    ``report`` (rows such as ``Mapped (%)``, ``Avg. coverage depth``, ``# SNPs`` ...).
    """
    out = _output_directory(output_directory)
    arguments, files = _common(contigs, threads, min_contig, False, False, None, False)
    if reference is not None:
        arguments += ["-r", str(_cli.existing_file(reference, "reference"))]
    if (pe1 is None) != (pe2 is None):
        raise ConfigurationError("pe1 and pe2 must be given together.")
    options = (
        ("--pe1", pe1),
        ("--pe2", pe2),
        ("--single", single),
        ("--pacbio", pacbio),
        ("--nanopore", nanopore),
        ("--ref-bam", reference_bam),
    )
    given = False
    for flag, value in options:
        if value is not None:
            arguments += [flag, str(_cli.existing_file(value, flag.lstrip("-")))]
            given = True
    if contig_bam_files is not None:
        bams = [str(path) for path in _cli.existing_files(contig_bam_files, "contig_bam_files")]
        if len(bams) != len(files):
            raise ConfigurationError("contig_bam_files needs one BAM per contig file.")
        arguments += ["--bam", ",".join(bams)]
        given = True
    if not given:
        raise ConfigurationError(
            "give reads (pe1/pe2, single, pacbio, nanopore) or alignments (reference_bam, contig_bam_files)."
        )
    if not _cli.flag(detect_variants, "detect_variants"):
        arguments += ["--no-snps", "--no-sv"]
    return _run(_TOOL, arguments, files, out, quast_path, timeout_seconds, "map_reads_to_assembly")


def metaquast_assess_metagenome(
    contigs: Sequence[str],
    output_directory: str,
    *,
    references: Sequence[str] | None = None,
    references_list: str | None = None,
    max_references: int = 50,
    use_blast_search: bool = True,
    unique_mapping: bool = False,
    min_contig: int = 500,
    threads: int = 1,
    quast_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate metagenome assemblies against several reference genomes, optionally found by BLAST (``metaquast.py``).

    ``references`` are reference genome FASTA files (``-r``, comma-joined); MetaQUAST aligns the contigs to all of them and
    reports per-reference statistics. Without references, or in addition to them, it searches the SILVA 16S database and
    downloads matching NCBI genomes (``use_blast_search=True``, which needs network access; ``False`` is ``--max-ref-number 0``
    when no references are given). ``references_list`` is a text file of reference names to download, ``max_references`` caps the
    number of downloaded genomes and ``unique_mapping`` (``--unique-mapping``) assigns every contig to one reference. The
    combined report is in ``report`` and the per-reference reports are in the ``runs_per_reference`` folder.
    """
    out = _output_directory(output_directory)
    arguments, files = _common(contigs, threads, min_contig, False, False, None, False)
    if references:
        paths = [str(path) for path in _cli.existing_files(references, "references")]
        arguments += ["-r", ",".join(paths)]
    if references_list is not None:
        arguments += [
            "--references-list",
            str(_cli.existing_file(references_list, "references_list")),
        ]
    if (
        not _cli.flag(use_blast_search, "use_blast_search")
        and not references
        and references_list is None
    ):
        raise ConfigurationError("without BLAST search give references or references_list.")
    if not use_blast_search:
        arguments += ["--max-ref-number", "0"]
    else:
        arguments += [
            "--max-ref-number",
            str(_cli.integer(max_references, "max_references", minimum=0)),
        ]
    if _cli.flag(unique_mapping, "unique_mapping"):
        arguments.append("--unique-mapping")
    return _run(
        _META_TOOL,
        arguments,
        files,
        out,
        quast_path,
        timeout_seconds,
        "metagenome_assembly_evaluation",
    )


__all__ = [
    "metaquast_assess_metagenome",
    "quast_assess_assembly",
    "quast_map_reads",
    "quast_predict_genes",
]
_register_public_verbose(globals())
