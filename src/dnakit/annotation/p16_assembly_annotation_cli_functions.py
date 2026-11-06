"""Assembly, annotation lift-over and protein alignment facades (pack P16): Flye, Liftoff, miniprot.

``assemble_long_reads_flye`` (including metaFlye and scaffolding), ``polish_assembly_flye``, ``lift_over_annotation_liftoff``
(extra copies, polishing), ``build_miniprot_index`` and ``align_proteins_to_genome_miniprot`` (PAF, GFF3 or translated
output). Programs run without a shell through ``dnakit.bioinformatics._p16_cli``; results are JSON-serialisable
dictionaries with the produced files and a bounded log tail.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli

_FLYE = cli.Tool(
    "flye",
    "flye",
    "Install Flye (conda install -c bioconda flye) and put it on PATH or pass flye_path.",
    "2.9.5",
    "https://doi.org/10.1038/s41587-019-0072-8",
    ("--version",),
    None,
)
_LIFTOFF = cli.Tool(
    "liftoff",
    "liftoff",
    "Install Liftoff (conda install -c bioconda liftoff; needs minimap2) and put it on PATH or pass liftoff_path.",
    "1.6.3",
    "https://doi.org/10.1093/bioinformatics/btaa1016",
    ("--version",),
    None,
)
_MINIPROT = cli.Tool(
    "miniprot",
    "miniprot",
    "Install miniprot (conda install -c bioconda miniprot) and put it on PATH or pass miniprot_path.",
    "0.18",
    "https://doi.org/10.1093/bioinformatics/btad014",
    ("--version",),
    None,
)
_READ_TYPES = ("pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq")


def assemble_long_reads_flye(
    reads: str,
    output_directory: str,
    *,
    read_type: str = "nano-raw",
    genome_size: str | None = None,
    metagenome: bool = False,
    scaffold: bool = False,
    iterations: int = 1,
    min_overlap: int | None = None,
    threads: int = 1,
    flye_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """De novo assembly of long reads with Flye, optionally in metagenome mode and with scaffolding.

    ``read_type`` is one of pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, nano-hq. ``metagenome=True``
    (``--meta``) is metaFlye for uneven coverage; ``scaffold=True`` (``--scaffold``) joins contigs into gapped
    scaffolds along the repeat graph. ``genome_size`` (``--genome-size``, for example ``"5m"``) is optional in recent
    Flye. ``iterations`` is the number of polishing rounds (``--iterations``). ``output_directory`` must be new; the
    assembly is ``assembly.fasta`` with ``assembly_info.txt`` (parsed into ``contigs``).
    """
    source = cli.existing_file(reads, "reads")
    out_dir = cli.new_directory(output_directory, "output_directory")
    arguments = [
        f"--{cli.choice(read_type, _READ_TYPES, 'read_type')}",
        str(source),
        "--out-dir",
        str(out_dir),
        "--threads",
        str(cli.threads_value(threads)),
        "--iterations",
        str(cli.integer(iterations, "iterations", minimum=0)),
    ]
    if genome_size is not None:
        arguments += ["--genome-size", cli.text(genome_size, "genome_size")]
    if cli.flag(metagenome, "metagenome"):
        arguments.append("--meta")
    if cli.flag(scaffold, "scaffold"):
        arguments.append("--scaffold")
    if min_overlap is not None:
        arguments += ["--min-overlap", str(cli.integer(min_overlap, "min_overlap", minimum=1))]
    assembly, info = out_dir / "assembly.fasta", out_dir / "assembly_info.txt"
    record = cli.invoke(
        _FLYE, arguments, cwd=out_dir, expected=[assembly, info], explicit=flye_path
    )
    contigs = cli.read_table(
        info, comment="#", limit=1000, text_columns=frozenset({"seq_name", "circ.", "graph_path"})
    )
    return cli.finish(
        record,
        "assemble_long_reads_flye",
        metagenome=metagenome,
        scaffold=scaffold,
        contigs=contigs,
    )


def polish_assembly_flye(
    draft: str,
    reads: str,
    output_directory: str,
    *,
    read_type: str = "nano-raw",
    iterations: int = 1,
    threads: int = 1,
    flye_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Polish a given draft assembly with long reads (standalone Flye polisher, ``flye --polish-target``).

    ``reads`` may be FASTA/FASTQ or an aligned BAM; the polished contigs are ``polished_<iterations>.fasta``.
    """
    contigs = cli.existing_file(draft, "draft")
    source = cli.existing_file(reads, "reads")
    out_dir = cli.new_directory(output_directory, "output_directory")
    rounds = cli.integer(iterations, "iterations", minimum=1)
    arguments = [
        "--polish-target",
        str(contigs),
        f"--{cli.choice(read_type, _READ_TYPES, 'read_type')}",
        str(source),
        "--out-dir",
        str(out_dir),
        "--iterations",
        str(rounds),
        "--threads",
        str(cli.threads_value(threads)),
    ]
    polished = out_dir / f"polished_{rounds}.fasta"
    record = cli.invoke(_FLYE, arguments, cwd=out_dir, expected=[polished], explicit=flye_path)
    return cli.finish(record, "polish_assembly_flye", polished=str(polished))


def lift_over_annotation_liftoff(
    target_fasta: str,
    reference_fasta: str,
    reference_annotation: str,
    output_gff: str,
    *,
    unmapped_path: str | None = None,
    extra_copies: bool = False,
    polish: bool = False,
    min_coverage: float = 0.5,
    min_identity: float = 0.5,
    copy_identity: float = 0.95,
    feature_types: str | None = None,
    threads: int = 1,
    liftoff_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift a GFF/GTF gene annotation from a reference to a target assembly (``liftoff``).

    Genes are aligned with minimap2 and mapped exon by exon. ``extra_copies=True`` (``-copies``) searches for additional
    gene copies with sequence identity at least ``copy_identity`` (``-sc``); ``polish=True`` (``-polish``) re-aligns exons to
    restore CDS start and stop codons. ``min_coverage`` (``-a``) and ``min_identity`` (``-s``) are Liftoff's thresholds;
    ``feature_types`` is a feature-type list file (``-f``). Unmapped gene ids go to ``unmapped_path`` (default
    ``<output_gff>.unmapped.txt``).
    """
    target = cli.existing_file(target_fasta, "target_fasta")
    reference = cli.existing_file(reference_fasta, "reference_fasta")
    annotation = cli.existing_file(reference_annotation, "reference_annotation")
    out = cli.new_file(output_gff, "output_gff")
    unmapped = (
        cli.new_file(unmapped_path, "unmapped_path")
        if unmapped_path
        else out.with_name(out.name + ".unmapped.txt")
    )
    arguments = [
        "-g",
        str(annotation),
        "-o",
        str(out),
        "-u",
        str(unmapped),
        "-a",
        str(cli.number(min_coverage, "min_coverage", minimum=0.0, maximum=1.0)),
        "-s",
        str(cli.number(min_identity, "min_identity", minimum=0.0, maximum=1.0)),
        "-p",
        str(cli.threads_value(threads)),
        "-dir",
        str(out.with_name(out.name + ".intermediate")),
    ]
    if cli.flag(extra_copies, "extra_copies"):
        arguments += [
            "-copies",
            "-sc",
            str(cli.number(copy_identity, "copy_identity", minimum=0.0, maximum=1.0)),
        ]
    if cli.flag(polish, "polish"):
        arguments.append("-polish")
    if feature_types is not None:
        arguments += ["-f", str(cli.existing_file(feature_types, "feature_types"))]
    arguments += [str(target), str(reference)]
    record = cli.invoke(
        _LIFTOFF, arguments, cwd=out.parent, expected=[out, unmapped], explicit=liftoff_path
    )
    genes = sum(
        1
        for line in cli.open_text(out)
        if not line.startswith("#") and line.split("\t")[2:3] == ["gene"]
    )
    return cli.finish(
        record,
        "lift_over_annotation_liftoff",
        lifted_genes=genes,
        unmapped_features=cli.count_lines(unmapped),
        extra_copies=extra_copies,
        polish=polish,
    )


def build_miniprot_index(
    genome: str,
    index_path: str,
    *,
    kmer: int = 6,
    translation_table: int = 1,
    threads: int = 1,
    miniprot_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a reusable miniprot genome index (``miniprot -d FILE ref.fna``)."""
    ref = cli.existing_file(genome, "genome")
    index = cli.new_file(index_path, "index_path")
    arguments = [
        "-t",
        str(cli.threads_value(threads)),
        "-k",
        str(cli.integer(kmer, "kmer", minimum=1, maximum=10)),
        "-T",
        str(cli.integer(translation_table, "translation_table", minimum=1, maximum=33)),
        "-d",
        str(index),
        str(ref),
    ]
    record = cli.invoke(
        _MINIPROT, arguments, cwd=index.parent, expected=[index], explicit=miniprot_path
    )
    return cli.finish(record, "build_miniprot_index", index=str(index))


def align_proteins_to_genome_miniprot(
    genome_or_index: str,
    proteins: str,
    output_path: str,
    *,
    output_format: str = "paf",
    translate: bool = False,
    include_alignment: bool = False,
    no_splicing: bool = False,
    max_intron_size: int | None = None,
    splice_model: int | None = None,
    translation_table: int = 1,
    threads: int = 1,
    miniprot_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align protein sequences to a genome, spliced and frameshift-aware (``miniprot``).

    ``genome_or_index`` is a genome FASTA or an index from :func:`build_miniprot_index`. ``output_format`` is ``paf``
    (default), ``gff`` (GFF3, ``--gff``) or ``gtf`` (``--gtf``). ``translate`` adds the translated protein of each
    alignment (``--trans``, skipping frameshifts) and ``include_alignment`` the residue alignment (``--aln``); both apply to
    the PAF/GFF3 output. ``no_splicing`` is ``-S``; ``max_intron_size`` is ``-G``; ``splice_model`` is ``-j``
    (2 vertebrate/insect, 1 general, 0 none). The result file is written to ``output_path``.
    """
    ref = cli.existing_file(genome_or_index, "genome_or_index")
    prot = cli.existing_file(proteins, "proteins")
    target = cli.new_file(output_path, "output_path")
    kind = cli.choice(output_format, ("paf", "gff", "gtf"), "output_format")
    arguments = ["-t", str(cli.threads_value(threads))]
    if kind == "gff":
        arguments.append("--gff")
    elif kind == "gtf":
        arguments.append("--gtf")
    if cli.flag(translate, "translate"):
        arguments.append("--trans")
    if cli.flag(include_alignment, "include_alignment"):
        arguments.append("--aln")
    if cli.flag(no_splicing, "no_splicing"):
        arguments.append("-S")
    if max_intron_size is not None:
        arguments += ["-G", str(cli.integer(max_intron_size, "max_intron_size", minimum=1))]
    if splice_model is not None:
        arguments += ["-j", str(cli.integer(splice_model, "splice_model", minimum=0, maximum=2))]
    if ref.suffix != ".mpi":
        arguments += [
            "-T",
            str(cli.integer(translation_table, "translation_table", minimum=1, maximum=33)),
        ]
    arguments += [str(ref), str(prot)]
    record = cli.invoke(
        _MINIPROT, arguments, cwd=target.parent, stdout_to=target, explicit=miniprot_path
    )
    lines = [line for line in cli.open_text(target) if line.strip() and not line.startswith("##")]
    return cli.finish(
        record,
        "align_proteins_to_genome_miniprot",
        output_format=kind,
        lines=len(lines),
        head=[line.rstrip("\n") for line in lines[:20]],
    )


__all__ = [
    "align_proteins_to_genome_miniprot",
    "assemble_long_reads_flye",
    "build_miniprot_index",
    "lift_over_annotation_liftoff",
    "polish_assembly_flye",
]

_register_public_verbose(globals(), tuple(__all__))
del Path
