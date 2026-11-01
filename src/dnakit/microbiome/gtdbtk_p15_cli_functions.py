"""GTDB-Tk command-line facades: marker identification, alignment, placement, de novo tree, ANI and MSA tools.

GTDB-Tk (2.x, tested against the 2.8.0 command-line interface) classifies bacterial and archaeal genomes against the
Genome Taxonomy Database. The user installs GTDB-Tk and its programs (``conda install -c bioconda gtdbtk``) and downloads the
GTDB reference package, then passes its folder as ``database`` (DNAKit sets ``GTDBTK_DATA_PATH`` for the one process it
starts) or exports ``GTDBTK_DATA_PATH`` beforehand. Nothing is downloaded and no shell is used. Every function writes into a new
output folder and returns a JSON-serialisable dictionary with the status, backend, version, command, the files produced and
the parsed summary tables.

The end-to-end ``classify_wf`` command is ``dnakit.microbiome.external.classify_gtdbtk``; this module holds the other
GTDB-Tk commands.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="gtdbtk",
    executable="gtdbtk",
    hint=(
        "Install GTDB-Tk 2.4 or later (for example `conda install -c bioconda gtdbtk`) and download the GTDB reference package "
        "with `download-db.sh`; pass its folder as database or set GTDBTK_DATA_PATH. Put gtdbtk on PATH or pass gtdbtk_path."
    ),
    tested_version="2.8.0",
    citation="https://doi.org/10.1093/bioinformatics/btac672",
    version_arguments=("--version",),
    version_pattern=r"v?(\d+\.\d+\.\d+)",
)
_PROT_MODELS = ("JTT", "WAG", "LG")


def _environment(database: str | None) -> dict[str, str]:
    if database is not None:
        return {"GTDBTK_DATA_PATH": str(_cli.existing_directory(database, "database"))}
    if not os.environ.get("GTDBTK_DATA_PATH"):
        raise ConfigurationError("give database (the GTDB reference package folder) or set the GTDBTK_DATA_PATH variable.")
    return {}


def _genomes(genome_dir: str | None, batchfile: str | None) -> list[str]:
    if (genome_dir is None) == (batchfile is None):
        raise ConfigurationError("give exactly one of genome_dir and batchfile.")
    if genome_dir is not None:
        return ["--genome_dir", str(_cli.existing_directory(genome_dir, "genome_dir"))]
    return ["--batchfile", str(_cli.existing_file(batchfile, "batchfile"))]  # type: ignore[arg-type]


def _common(prefix: str, cpus: int, extension: str | None = None, debug: bool = False) -> list[str]:
    arguments = ["--prefix", _cli.text(prefix, "prefix"), "--cpus", str(_cli.integer(cpus, "cpus", minimum=1, maximum=512))]
    if extension is not None:
        arguments += ["--extension", _cli.text(extension, "extension")]
    if _cli.flag(debug, "debug"):
        arguments.append("--debug")
    return arguments


def _tables(directory: Path, pattern: str, limit: int = 2000) -> dict[str, list[dict[str, Any]]]:
    return {
        str(path.relative_to(directory)): _cli.read_table(path, limit=limit)
        for path in sorted(directory.rglob(pattern))
        if path.is_file() and path.stat().st_size > 0
    }


def _files(directory: Path, patterns: Sequence[str]) -> dict[str, list[str]]:
    return {pattern: sorted(str(path) for path in directory.rglob(pattern) if path.is_file()) for pattern in patterns}


def _run(
    arguments: Sequence[str], out: Path, database: str | None, gtdbtk_path: str | None, timeout_seconds: float,
    expected: Sequence[Path] = (),
) -> dict[str, Any]:
    return _cli.invoke(
        _TOOL, arguments, cwd=out.parent if out.is_dir() and not out.samefile(out.parent) else out, expected=expected,
        explicit=gtdbtk_path, timeout_seconds=timeout_seconds, environment=_environment(database),
    )


def gtdbtk_identify_markers(
    out_dir: str,
    *,
    genome_dir: str | None = None,
    batchfile: str | None = None,
    database: str | None = None,
    extension: str = "fna",
    prefix: str = "gtdbtk",
    genes: bool = False,
    write_single_copy_genes: bool = False,
    cpus: int = 1,
    force: bool = False,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Identify the GTDB archaeal (53) and bacterial (120) marker genes in genome FASTA files (``gtdbtk identify``).

    Genes are called with Prodigal (``genes=True`` treats the inputs as protein FASTA of already predicted genes) and the
    markers are found with HMMER. Give either ``genome_dir`` (all files with ``extension``) or ``batchfile`` (tab-separated
    path and genome name). ``write_single_copy_genes`` also writes the single-copy marker sequences; ``force`` carries on when
    a genome fails. Output goes to the new folder ``out_dir``; the result lists the produced files and the parsed
    ``*.markers_summary.tsv`` tables (``markers_summary``) and the translation table summary.
    """
    out = _cli.new_directory(out_dir, "out_dir")
    arguments = ["identify", *_genomes(genome_dir, batchfile), "--out_dir", str(out), *_common(prefix, cpus, extension)]
    for enabled, option in ((genes, "--genes"), (write_single_copy_genes, "--write_single_copy_genes"), (force, "--force")):
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _run(arguments, out, database, gtdbtk_path, timeout_seconds)
    return _cli.finish(
        record, "identify_gtdb_marker_genes", out_dir=str(out),
        markers_summary=_tables(out, "*.markers_summary.tsv"), translation_table_summary=_tables(out, "*.translation_table_summary.tsv"),
        failed_genomes=_tables(out, "*.failed_genomes.tsv"), files=_files(out / "identify", ["*.tsv"]),
    )


def gtdbtk_align_markers(
    out_dir: str,
    identify_dir: str,
    *,
    database: str | None = None,
    prefix: str = "gtdbtk",
    skip_gtdb_refs: bool = False,
    taxa_filter: str | None = None,
    min_perc_aa: float = 10.0,
    custom_msa_filters: bool = False,
    skip_trimming: bool = False,
    cols_per_gene: int = 42,
    min_consensus: float = 25.0,
    max_consensus: float = 95.0,
    min_perc_taxa: float = 50.0,
    rnd_seed: int | None = None,
    cpus: int = 1,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align the identified marker genes of the genomes into a masked multiple sequence alignment (``gtdbtk align``).

    ``identify_dir`` is the output folder of ``gtdbtk_identify_markers``. The user genomes are aligned with the GTDB reference
    markers (``skip_gtdb_refs`` leaves the references out; ``taxa_filter`` keeps only references of the named taxa) and the
    alignment is trimmed with the GTDB mask unless ``skip_trimming``; ``custom_msa_filters`` instead applies
    ``cols_per_gene``, ``min_consensus``, ``max_consensus`` and ``min_perc_taxa``. ``min_perc_aa`` drops genomes with too few
    aligned amino acids. The alignments (``*.msa.fasta`` and ``*.user_msa.fasta``, possibly gzip-compressed) are listed in
    ``files``; ``filtered_genomes`` holds the genomes removed by the filters.
    """
    out = _cli.new_directory(out_dir, "out_dir")
    identify = _cli.existing_directory(identify_dir, "identify_dir")
    if custom_msa_filters and skip_trimming:
        raise ConfigurationError("custom_msa_filters and skip_trimming cannot be used together.")
    arguments = ["align", "--identify_dir", str(identify), "--out_dir", str(out), *_common(prefix, cpus)]
    arguments += ["--min_perc_aa", str(_cli.number(min_perc_aa, "min_perc_aa", minimum=0.0, maximum=100.0))]
    arguments += ["--cols_per_gene", str(_cli.integer(cols_per_gene, "cols_per_gene", minimum=1))]
    arguments += ["--min_consensus", str(_cli.number(min_consensus, "min_consensus", minimum=0.0, maximum=100.0))]
    arguments += ["--max_consensus", str(_cli.number(max_consensus, "max_consensus", minimum=0.0, maximum=100.0))]
    arguments += ["--min_perc_taxa", str(_cli.number(min_perc_taxa, "min_perc_taxa", minimum=0.0, maximum=100.0))]
    if taxa_filter is not None:
        arguments += ["--taxa_filter", _cli.text(taxa_filter, "taxa_filter")]
    if rnd_seed is not None:
        arguments += ["--rnd_seed", str(_cli.integer(rnd_seed, "rnd_seed"))]
    for enabled, option in ((skip_gtdb_refs, "--skip_gtdb_refs"), (custom_msa_filters, "--custom_msa_filters"), (skip_trimming, "--skip_trimming")):
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _run(arguments, out, database, gtdbtk_path, timeout_seconds)
    return _cli.finish(
        record, "align_gtdb_marker_genes", out_dir=str(out),
        files=_files(out / "align", ["*.msa.fasta*", "*.user_msa.fasta*", "*.filtered.tsv"]), filtered_genomes=_tables(out, "*.filtered.tsv"),
    )


def gtdbtk_classify_placement(
    out_dir: str,
    align_dir: str,
    *,
    genome_dir: str | None = None,
    batchfile: str | None = None,
    database: str | None = None,
    extension: str = "fna",
    prefix: str = "gtdbtk",
    full_tree: bool = False,
    place_species: bool = False,
    min_af: float | None = None,
    skip_pplacer: bool = False,
    pplacer_cpus: int | None = None,
    cpus: int = 1,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify genomes by pplacer placement of a precomputed marker alignment (``gtdbtk classify``).

    ``align_dir`` is the output folder of ``gtdbtk_align_markers``; ``genome_dir`` or ``batchfile`` names the same genomes
    again (the ANI step needs the sequences). Genomes are placed in the GTDB reference tree (``full_tree`` uses the unsplit
    bacterial tree, which needs much more memory), the taxonomy is assigned from the placement, relative evolutionary
    divergence and ANI to the closest representative (``min_af`` sets the minimum alignment fraction, ``place_species`` runs the
    skani species screening first, ``skip_pplacer`` stops after it). Result keys: ``summary`` (the parsed
    ``*.summary.tsv`` tables with user genome, classification, closest reference, ANI and warnings) and ``trees``.
    """
    out = _cli.new_directory(out_dir, "out_dir")
    align = _cli.existing_directory(align_dir, "align_dir")
    arguments = ["classify", *_genomes(genome_dir, batchfile), "--align_dir", str(align), "--out_dir", str(out), *_common(prefix, cpus, extension)]
    if pplacer_cpus is not None:
        arguments += ["--pplacer_cpus", str(_cli.integer(pplacer_cpus, "pplacer_cpus", minimum=1))]
    if min_af is not None:
        arguments += ["--min_af", str(_cli.number(min_af, "min_af", minimum=0.0, maximum=1.0))]
    for enabled, option in ((full_tree, "--full_tree"), (place_species, "--place_species"), (skip_pplacer, "--skip_pplacer")):
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _run(arguments, out, database, gtdbtk_path, timeout_seconds)
    return _cli.finish(
        record, "classify_gtdb_by_placement", out_dir=str(out), summary=_tables(out, "*.summary.tsv"),
        trees=_files(out / "classify", ["*.classify.tree*"]),
    )


def gtdbtk_de_novo_workflow(
    out_dir: str,
    outgroup_taxon: str,
    *,
    domain: str = "bacteria",
    genome_dir: str | None = None,
    batchfile: str | None = None,
    database: str | None = None,
    extension: str = "fna",
    prefix: str = "gtdbtk",
    skip_gtdb_refs: bool = False,
    taxa_filter: str | None = None,
    min_perc_aa: float = 10.0,
    custom_msa_filters: bool = False,
    cols_per_gene: int = 42,
    min_consensus: float = 25.0,
    max_consensus: float = 95.0,
    min_perc_taxa: float = 50.0,
    prot_model: str = "WAG",
    no_support: bool = False,
    gamma: bool = False,
    rnd_seed: int | None = None,
    genes: bool = False,
    write_single_copy_genes: bool = False,
    keep_intermediates: bool = False,
    force: bool = False,
    cpus: int = 1,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer a de novo genome tree rooted and decorated with GTDB taxonomy (``gtdbtk de_novo_wf``).

    Runs identify, align, infer (FastTree with ``prot_model`` ``JTT``, ``WAG`` or ``LG``; ``no_support`` skips support values,
    ``gamma`` adds gamma rate heterogeneity), root (with the GTDB taxon ``outgroup_taxon`` such as ``p__Patescibacteria``)
    and decorate. ``domain`` is ``bacteria`` (``--bacteria``) or ``archaea`` (``--archaea``). The MSA options are those of
    ``gtdbtk_align_markers``. The Newick trees (``*.decorated.tree``, ``*.rooted.tree``) are listed in ``trees`` together with
    the files of every stage.
    """
    out = _cli.new_directory(out_dir, "out_dir")
    kingdom = {"bacteria": "--bacteria", "archaea": "--archaea"}[_cli.choice(domain, ("bacteria", "archaea"), "domain")]
    arguments = ["de_novo_wf", *_genomes(genome_dir, batchfile), kingdom, "--outgroup_taxon", _cli.text(outgroup_taxon, "outgroup_taxon")]
    arguments += ["--out_dir", str(out), *_common(prefix, cpus, extension)]
    arguments += ["--min_perc_aa", str(_cli.number(min_perc_aa, "min_perc_aa", minimum=0.0, maximum=100.0))]
    arguments += ["--cols_per_gene", str(_cli.integer(cols_per_gene, "cols_per_gene", minimum=1))]
    arguments += ["--min_consensus", str(_cli.number(min_consensus, "min_consensus", minimum=0.0, maximum=100.0))]
    arguments += ["--max_consensus", str(_cli.number(max_consensus, "max_consensus", minimum=0.0, maximum=100.0))]
    arguments += ["--min_perc_taxa", str(_cli.number(min_perc_taxa, "min_perc_taxa", minimum=0.0, maximum=100.0))]
    arguments += ["--prot_model", _cli.choice(prot_model, _PROT_MODELS, "prot_model")]
    if taxa_filter is not None:
        arguments += ["--taxa_filter", _cli.text(taxa_filter, "taxa_filter")]
    if rnd_seed is not None:
        arguments += ["--rnd_seed", str(_cli.integer(rnd_seed, "rnd_seed"))]
    flags = (
        (skip_gtdb_refs, "--skip_gtdb_refs"), (custom_msa_filters, "--custom_msa_filters"), (no_support, "--no_support"),
        (gamma, "--gamma"), (genes, "--genes"), (write_single_copy_genes, "--write_single_copy_genes"),
        (keep_intermediates, "--keep_intermediates"), (force, "--force"),
    )
    for enabled, option in flags:
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _run(arguments, out, database, gtdbtk_path, timeout_seconds)
    return _cli.finish(
        record, "gtdb_de_novo_workflow", out_dir=str(out), trees=_files(out, ["*.decorated.tree", "*.rooted.tree", "*.unrooted.tree"]),
        files=_files(out, ["*.tsv", "*.msa.fasta*"]),
    )


def gtdbtk_ani_to_representatives(
    out_dir: str,
    *,
    genome_dir: str | None = None,
    batchfile: str | None = None,
    database: str | None = None,
    extension: str = "fna",
    prefix: str = "gtdbtk",
    min_af: float | None = None,
    cpus: int = 1,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average nucleotide identity of genomes to the GTDB representative genomes with skani (``gtdbtk ani_rep``).

    Each genome is compared with all GTDB species representatives; ``min_af`` is the minimum alignment fraction of a hit
    (GTDB-Tk's default when omitted). The parsed ``<prefix>.ani_summary.tsv`` rows (``ani_summary``) give for every genome its
    closest representative, ANI, alignment fraction and the GTDB taxonomy of that representative.
    """
    out = _cli.new_directory(out_dir, "out_dir")
    arguments = ["ani_rep", *_genomes(genome_dir, batchfile), "--out_dir", str(out), *_common(prefix, cpus, extension)]
    if min_af is not None:
        arguments += ["--min_af", str(_cli.number(min_af, "min_af", minimum=0.0, maximum=1.0))]
    summary = out / f"{prefix}.ani_summary.tsv"
    record = _run(arguments, out, database, gtdbtk_path, timeout_seconds, expected=[summary])
    return _cli.finish(record, "gtdb_ani_to_representatives", out_dir=str(out), ani_summary_file=str(summary), ani_summary=_cli.read_table(summary))


def gtdbtk_trim_msa(
    untrimmed_msa: str,
    output_path: str,
    *,
    mask_file: str | None = None,
    reference_mask: str | None = None,
    database: str | None = None,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim an untrimmed marker multiple sequence alignment with a column mask (``gtdbtk trim_msa``).

    Columns where the mask has ``1`` are kept. Use ``mask_file`` (a one-line text file of 0 and 1) or ``reference_mask``
    (``'bac'`` or ``'arc'``, the mask of the GTDB reference package, which needs ``database`` or ``GTDBTK_DATA_PATH``). The
    trimmed FASTA alignment is written to the new file ``output_path``.
    """
    source = _cli.existing_file(untrimmed_msa, "untrimmed_msa")
    target = _cli.new_file(output_path, "output_path")
    if (mask_file is None) == (reference_mask is None):
        raise ConfigurationError("give exactly one of mask_file and reference_mask.")
    arguments = ["trim_msa", "--untrimmed_msa", str(source), "--output", str(target)]
    if mask_file is not None:
        arguments += ["--mask_file", str(_cli.existing_file(mask_file, "mask_file"))]
    else:
        arguments += ["--reference_mask", _cli.choice(reference_mask, ("bac", "arc"), "reference_mask")]
    record = _cli.invoke(
        _TOOL, arguments, cwd=target.parent, expected=[target], explicit=gtdbtk_path, timeout_seconds=timeout_seconds,
        environment=_environment(database),
    )
    return _cli.finish(record, "trim_gtdb_msa", output=str(target), sequence_count=_cli.count_fastx_records(target))


def gtdbtk_export_reference_msa(
    output_path: str,
    *,
    domain: str = "bac",
    database: str | None = None,
    gtdbtk_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export the untrimmed GTDB reference marker alignment of the archaeal or bacterial domain (``gtdbtk export_msa``).

    ``domain`` is ``'bac'`` (bacterial, 120 markers) or ``'arc'`` (archaeal, 53 markers). The alignment of all GTDB reference
    genomes is written to the new file ``output_path``; the reference package comes from ``database`` or ``GTDBTK_DATA_PATH``.
    """
    target = _cli.new_file(output_path, "output_path")
    arguments = ["export_msa", "--domain", _cli.choice(domain, ("bac", "arc"), "domain"), "--output", str(target)]
    record = _cli.invoke(
        _TOOL, arguments, cwd=target.parent, expected=[target], explicit=gtdbtk_path, timeout_seconds=timeout_seconds,
        environment=_environment(database),
    )
    return _cli.finish(record, "export_gtdb_reference_msa", output=str(target), domain=domain, sequence_count=_cli.count_fastx_records(target))


__all__ = [
    "gtdbtk_align_markers",
    "gtdbtk_ani_to_representatives",
    "gtdbtk_classify_placement",
    "gtdbtk_de_novo_workflow",
    "gtdbtk_export_reference_msa",
    "gtdbtk_identify_markers",
    "gtdbtk_trim_msa",
]
_register_public_verbose(globals())
