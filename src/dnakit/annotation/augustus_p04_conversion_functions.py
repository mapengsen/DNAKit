"""AUGUSTUS gene-model, sequence and masking conversions (integration pack P04, part 3).

Each public function runs one AUGUSTUS 3.5.0 helper script (``gtf2gff.pl``, ``gff2gbSmallDNA.pl``,
``getAnnoFasta.pl``, ``filterGenes*.pl``, ``compare_masking.pl`` and others) with an explicit argument list and no
shell. The user installs AUGUSTUS (``conda install -c bioconda augustus``) and puts its ``bin`` folder on PATH or
passes ``augustus_bin_dir``. Outputs are new explicit paths and nothing is overwritten. Results are JSON-serialisable
dictionaries with status, backend, version, command, outputs, small parsed counts and a bounded log tail.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation import _augustus_p04_cli as _a
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError


def _in(path: object, name: str) -> Path:
    return _cli.existing_file(str(path), name)


def _out(path: object, name: str = "output_path") -> Path:
    return _cli.new_file(str(path), name)


def _stdin(
    operation: str,
    tool: str,
    input_path: str,
    output_path: str,
    arguments: list[str],
    directory: str | None,
    **extra: Any,
) -> dict[str, Any]:
    return _a.simple(
        operation,
        tool,
        arguments,
        directory=_a.bin_dir(directory),
        output=_out(output_path),
        stdin=_in(input_path, "input_path"),
        **extra,
    )


def augustus_add_names_to_gff3(
    gff3_path: str,
    output_path: str,
    *,
    only_gene_structure: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add ``Name=`` attributes to the gene and mRNA lines of an AUGUSTUS GFF3 (``add_name_to_gff3.pl``); ``only_gene_structure`` keeps gene, mRNA, CDS and exon only."""
    target = _out(output_path)
    arguments = [f"--out={target}"] + _a.switch(
        "--filter", _a.flag(only_gene_structure, "only_gene_structure")
    )
    return _a.simple(
        "add_names_to_gff3",
        "add_name_to_gff3.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        stdin=_in(gff3_path, "gff3_path"),
        to_stdout=False,
    )


def augustus_cegma_to_gff(
    cegma_gff_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a CEGMA output GFF to the CDS GTF used for AUGUSTUS training (``cegma2gff.pl``)."""
    return _a.simple(
        "cegma_to_gff",
        "cegma2gff.pl",
        [str(_in(cegma_gff_path, "cegma_gff_path"))],
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_scipio_gff_to_gff(
    scipio_gff_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a Scipio GFF to an AUGUSTUS GFF that ``augustus_gff_to_genbank`` can read (``scipiogff2gff.pl``)."""
    target = _out(output_path)
    return _a.simple(
        "scipio_gff_to_gff",
        "scipiogff2gff.pl",
        [f"--in={_in(scipio_gff_path, 'scipio_gff_path')}", f"--out={target}"],
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_scipio_yaml_to_gff(
    scipio_yaml_path: str,
    output_path: str,
    *,
    filter_status: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert Scipio YAML alignments to GFF (``yaml2gff.1.4.pl``); ``filter_status`` drops alignments with that status (for example ``incomplete``)."""
    arguments = (
        []
        if filter_status is None
        else [f"--filterstatus={_cli.text(filter_status, 'filter_status')}"]
    )
    return _stdin(
        "scipio_yaml_to_gff",
        "yaml2gff.1.4.pl",
        scipio_yaml_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_partition_gene_models(
    genome_path: str,
    gtf_path: str,
    output_path: str,
    *,
    min_cluster_size: int | None = None,
    max_flanking_region_size: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Partition single-exon CDS gene models into gene clusters and write them as a GenBank file, one cluster per LOCUS (``partition_gtf2gb.pl``)."""
    target = _out(output_path)
    arguments = [
        f"--genome={_in(genome_path, 'genome_path')}",
        f"--gtf={_in(gtf_path, 'gtf_path')}",
        f"--out={target}",
    ]
    arguments += _a.eq(
        "--minClusterSize", _a.optional_int(min_cluster_size, "min_cluster_size", minimum=1)
    )
    arguments += _a.eq(
        "--maxFlankingRegionSize",
        _a.optional_int(max_flanking_region_size, "max_flanking_region_size", minimum=1),
    )
    result = _a.simple(
        "partition_gene_models",
        "partition_gtf2gb.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )
    result["loci"] = sum(
        1
        for line in target.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("LOCUS")
    )
    return result


def augustus_cluster_gene_dependencies(
    dependency_path: str,
    clusters_path: str,
    *,
    threshold: int | None = None,
    max_cluster_size: int | None = None,
    transcript_format: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group genes into clusters without inter-cluster dependencies from a pair-dependency file (``clusterAndSplitGenes.pl``); clusters go to ``clusters_path``."""
    target = _out(clusters_path, "clusters_path")
    arguments = [f"--clusters={target}"] + _a.switch(
        "--txformat", _a.flag(transcript_format, "transcript_format")
    )
    arguments += _a.valued("-t", _a.optional_int(threshold, "threshold", minimum=1))
    arguments += _a.eq(
        "--maxclustersize", _a.optional_int(max_cluster_size, "max_cluster_size", minimum=1)
    )
    return _a.simple(
        "cluster_gene_dependencies",
        "clusterAndSplitGenes.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        stdin=_in(dependency_path, "dependency_path"),
        to_stdout=False,
    )


def augustus_compare_masking(
    softmasked_path: str,
    second_softmasked_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare the soft-masking of two versions of the same assembly (``compare_masking.pl``; both are held in memory). The report goes to ``output_path``."""
    arguments = [
        str(_in(softmasked_path, "softmasked_path")),
        str(_in(second_softmasked_path, "second_softmasked_path")),
    ]
    return _a.simple(
        "compare_masking",
        "compare_masking.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_merge_masking(
    softmasked_path: str,
    second_softmasked_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Soft-mask every base that is masked in either of two versions of the same assembly (``merge_masking.pl``)."""
    target = _out(output_path)
    arguments = [
        str(_in(softmasked_path, "softmasked_path")),
        str(_in(second_softmasked_path, "second_softmasked_path")),
        str(target),
    ]
    return _a.simple(
        "merge_masking",
        "merge_masking.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
        kind="fasta",
    )


def augustus_compute_flanking_region(
    gtf_path: str, output_path: str, *, augustus_bin_dir: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Average gene length and the flanking region derived from a GTF (``computeFlankingRegion.pl``); the script text is written to ``output_path`` and the numbers are parsed into ``numbers``."""
    result = _a.simple(
        "compute_flanking_region",
        "computeFlankingRegion.pl",
        [str(_in(gtf_path, "gtf_path"))],
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )
    import re

    text = Path(result["outputs"]["output"]).read_text(encoding="utf-8", errors="replace")
    result["numbers"] = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", text)]
    return result


_GENBANK_FILTERS = {
    "exclude_loci": "filterGenes.pl",
    "include_loci": "filterGenesIn.pl",
    "include_mrna_gtf": "filterGenesIn_mRNAname.pl",
    "exclude_mrna": "filterGenesOut_mRNAname.pl",
}


def augustus_filter_genbank_genes(
    genbank_path: str,
    names_path: str,
    output_path: str,
    *,
    mode: str = "exclude_loci",
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select gene models of a GenBank file by a gene-name list.

    ``mode`` ``exclude_loci`` drops the loci named in ``names_path`` (``filterGenes.pl``), ``include_loci`` keeps only
    them, first of identical loci (``filterGenesIn.pl``), ``include_mrna_gtf`` keeps the genes of a GTF (``names_path``
    is then a GTF, ``filterGenesIn_mRNAname.pl``) and ``exclude_mrna`` drops listed mRNA names (``filterGenesOut_mRNAname.pl``).
    """
    tool = _GENBANK_FILTERS[_cli.choice(mode, tuple(_GENBANK_FILTERS), "mode")]
    arguments = [str(_in(names_path, "names_path")), str(_in(genbank_path, "genbank_path"))]
    result = _a.simple(
        "filter_genbank_genes",
        tool,
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )
    result["loci"] = sum(
        1
        for line in Path(result["outputs"]["output"])
        .read_text(encoding="utf-8", errors="replace")
        .splitlines()
        if line.startswith("LOCUS")
    )
    return result


def augustus_check_gff_names_in_fasta(
    gff_path: str,
    genome_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that the sequence names of a GFF are FASTA headers of the genome and that it has 9 columns (``findGffNamesInFasta.pl``).

    ``names_missing`` is True when the GFF uses sequence names absent from the genome, ``wrong_column_count`` is True
    when a line does not have 9 columns (the script writes this flag to ``output_path``).
    """
    target = _out(output_path)
    arguments = [
        f"--gff={_in(gff_path, 'gff_path')}",
        f"--genome={_in(genome_path, 'genome_path')}",
        f"--out={target}",
    ]
    completed = _a.run(
        "findGffNamesInFasta.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        capture=True,
        expected=[target],
    )
    names_missing = _a.captured_text(completed).strip() == "1"
    wrong_columns = target.read_text(encoding="utf-8").strip() == "1"
    return _a.result(
        "check_gff_names_in_fasta",
        completed,
        outputs={"output": target},
        names_missing=names_missing,
        wrong_column_count=wrong_columns,
    )


def augustus_fix_joingenes_gtf(
    gtf_path: str, output_path: str, *, augustus_bin_dir: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Add the missing ``gene`` feature lines to a ``joingenes`` GTF (``fix_joingenes_gtf.pl``)."""
    return _stdin(
        "fix_joingenes_gtf", "fix_joingenes_gtf.pl", gtf_path, output_path, [], augustus_bin_dir
    )


def augustus_genbank_to_gff(
    genbank_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert GenBank training genes (with flanking-region LOCUS names such as ``chr1_100-900``) to genome-scale GFF (``gbSmallDNA2gff.pl``)."""
    return _a.simple(
        "genbank_to_gff",
        "gbSmallDNA2gff.pl",
        [str(_in(genbank_path, "genbank_path"))],
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_gbrowse_old_to_gff3(
    prediction_path: str,
    output_path: str,
    *,
    hints: bool = False,
    source: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert legacy GBrowse-format AUGUSTUS predictions to GFF3 (``gbrowseold2gff3.pl``); ``hints`` reads hints instead of genes."""
    arguments = _a.switch("--hints", _a.flag(hints, "hints")) + _a.eq(
        "--source", None if source is None else _cli.text(source, "source")
    )
    return _stdin(
        "gbrowse_old_to_gff3",
        "gbrowseold2gff3.pl",
        prediction_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_extract_annotation_fasta(
    prediction_path: str,
    output_dir: str,
    *,
    genome_path: str | None = None,
    chop_cds: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write protein and coding sequences from AUGUSTUS output comments (``getAnnoFasta.pl``).

    The prediction must contain ``# protein sequence`` and ``# coding sequence`` blocks (``protein_sequences=True`` and
    ``coding_sequences=True`` in ``augustus_predict_genes``). With ``genome_path`` the exon sequences (``.cdsexons``) and
    the mRNA including UTRs (``.mrna``) are written too. Files are ``augustus.aa``, ``.codingseq`` and so on in the new
    ``output_dir``.
    """
    source = _in(prediction_path, "prediction_path")
    out = _cli.new_directory(output_dir, "output_dir")
    arguments: list[str] = []
    if genome_path is not None:
        arguments.append(f"--seqfile={_in(genome_path, 'genome_path')}")
    arguments += _a.switch("--chop_cds", _a.flag(chop_cds, "chop_cds"))
    with tempfile.TemporaryDirectory(prefix="dnakit-augustus-fasta-") as scratch:
        copy = Path(scratch) / "augustus.gff"
        shutil.copyfile(source, copy)
        completed = _a.run(
            "getAnnoFasta.pl",
            [*arguments, str(copy)],
            directory=_a.bin_dir(augustus_bin_dir),
            cwd=Path(scratch),
            capture=True,
        )
        _a.captured_text(completed)
        produced = {}
        for suffix in ("aa", "codingseq", "cdsexons", "mrna"):
            path = Path(scratch) / f"augustus.{suffix}"
            if path.is_file():
                shutil.move(str(path), out / path.name)
                produced[suffix] = out / path.name
    counts = {suffix: _a.fasta_count(path) for suffix, path in produced.items()}
    return _a.result(
        "extract_annotation_fasta",
        completed,
        outputs={name: path for name, path in produced.items()},
        records=counts,
    )


def augustus_extract_cds_sequences(
    genome_path: str,
    annotation_path: str,
    output_stem: str,
    *,
    method: str = "getAnnoFastaFromJoingenes",
    annotation_format: str = "gtf",
    translation_table: int | None = None,
    drop_invalid_stops: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract coding sequences and translated proteins of CDS models in a GTF/GFF3 plus genome FASTA.

    ``method="getAnnoFastaFromJoingenes"`` (``getAnnoFastaFromJoingenes.py``) writes ``<stem>.codingseq`` and ``<stem>.aa``
    (GTF or GFF3, ``translation_table``, ``drop_invalid_stops``); ``"stringtie2fa"`` (``stringtie2fa.py``, GTF only) is the
    earlier variant with the same outputs; ``"gtf2aa"`` (``gtf2aa.pl``) writes only the protein FASTA ``<stem>.aa``.
    """
    genome = _in(genome_path, "genome_path")
    annotation = _in(annotation_path, "annotation_path")
    chosen = _cli.choice(method, ("getAnnoFastaFromJoingenes", "stringtie2fa", "gtf2aa"), "method")
    stem = Path(output_stem).expanduser().resolve()
    if not stem.parent.is_dir():
        raise ConfigurationError(f"The folder for output_stem does not exist: {stem.parent}")
    aa = stem.with_name(stem.name + ".aa")
    coding = stem.with_name(stem.name + ".codingseq")
    expected = [aa] if chosen == "gtf2aa" else [aa, coding]
    for path in expected:
        if path.exists():
            raise ConfigurationError(f"{path} already exists.")
    directory = _a.bin_dir(augustus_bin_dir)
    if chosen == "gtf2aa":
        tool, arguments = "gtf2aa.pl", [str(genome), str(annotation), str(aa)]
    else:
        fmt = _cli.choice(annotation_format, ("gtf", "gff3"), "annotation_format")
        tool = (
            "getAnnoFastaFromJoingenes.py"
            if chosen == "getAnnoFastaFromJoingenes"
            else "stringtie2fa.py"
        )
        if chosen == "stringtie2fa" and fmt != "gtf":
            raise ConfigurationError("stringtie2fa accepts GTF only.")
        arguments = [
            "-g",
            str(genome),
            "-o",
            str(stem),
            "-f" if fmt == "gtf" else "-3",
            str(annotation),
        ]
        if chosen == "getAnnoFastaFromJoingenes":
            if translation_table is not None:
                arguments += [
                    "-t",
                    str(
                        _cli.integer(translation_table, "translation_table", minimum=1, maximum=33)
                    ),
                ]
            if _a.flag(drop_invalid_stops, "drop_invalid_stops"):
                arguments += ["-s", "True"]
    completed = _a.run(
        tool, arguments, directory=directory, cwd=stem.parent, capture=True, expected=expected
    )
    _a.captured_text(completed)
    outputs = {"protein": aa} if chosen == "gtf2aa" else {"protein": aa, "coding": coding}
    return _a.result(
        "extract_cds_sequences",
        completed,
        outputs=outputs,
        method=chosen,
        records={name: _a.fasta_count(path) for name, path in outputs.items()},
    )


def augustus_extract_mrna_sequences(
    gff_path: str,
    genome_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Splice mRNA sequences of GFF gene models out of the genome (``gffGetmRNA.pl``) into the FASTA ``output_path``."""
    target = _out(output_path)
    arguments = [f"--genome={_in(genome_path, 'genome_path')}", f"--mrna={target}"]
    return _a.simple(
        "extract_mrna_sequences",
        "gffGetmRNA.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        stdin=_in(gff_path, "gff_path"),
        to_stdout=False,
        kind="fasta",
    )


def augustus_genbank_loci(
    genbank_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the tab-separated transcript-to-locus list of a GenBank training file (``get_loci_from_gb.pl``).

    The script prints the list to standard output and leaves its ``--out`` file empty (upstream bug); this facade
    points ``--out`` at a scratch file and writes the printed list to the new ``output_path``.
    """
    target = _out(output_path)
    with tempfile.TemporaryDirectory(prefix="dnakit-augustus-loci-") as scratch:
        arguments = [
            f"--gb={_in(genbank_path, 'genbank_path')}",
            f"--out={Path(scratch) / 'unused.lst'}",
        ]
        return _a.simple(
            "genbank_loci",
            "get_loci_from_gb.pl",
            arguments,
            directory=_a.bin_dir(augustus_bin_dir),
            output=target,
        )


def augustus_gff_to_genbank(
    gff_path: str,
    genome_path: str,
    flanking_dna: int,
    output_path: str,
    *,
    exclude_genes_path: str | None = None,
    include_genes_path: str | None = None,
    allow_overlaps: bool = False,
    connected: bool = False,
    keep_softmasking: bool = False,
    hardmask: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert GFF gene models and a genome FASTA to a minimal GenBank training file with flanking DNA (``gff2gbSmallDNA.pl``)."""
    target = _out(output_path)
    arguments = [
        str(_in(gff_path, "gff_path")),
        str(_in(genome_path, "genome_path")),
        str(_cli.integer(flanking_dna, "flanking_dna", minimum=0)),
        str(target),
    ]
    if exclude_genes_path is not None:
        arguments.append(f"--bad={_in(exclude_genes_path, 'exclude_genes_path')}")
    if include_genes_path is not None:
        arguments.append(f"--good={_in(include_genes_path, 'include_genes_path')}")
    arguments += _a.switch("--overlap", _a.flag(allow_overlaps, "allow_overlaps")) + _a.switch(
        "--connected", _a.flag(connected, "connected")
    )
    arguments += _a.switch(
        "--softmasked", _a.flag(keep_softmasking, "keep_softmasking")
    ) + _a.switch("--hardmask", _a.flag(hardmask, "hardmask"))
    result = _a.simple(
        "gff_to_genbank",
        "gff2gbSmallDNA.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )
    result["loci"] = sum(
        1
        for line in target.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("LOCUS")
    )
    return result


def augustus_convert_genepred(
    genepred_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert UCSC ``genePredExt`` rows to the bigGenePred-ready layout (``gp2othergp.pl``)."""
    return _stdin(
        "convert_genepred", "gp2othergp.pl", genepred_path, output_path, [], augustus_bin_dir
    )


def augustus_gtf_to_gff(
    gtf_path: str,
    output_path: str,
    *,
    gff3: bool = False,
    print_exons: bool = False,
    print_utrs: bool = False,
    print_introns: bool = False,
    include_stop_in_cds: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert AUGUSTUS GTF gene models to GFF (``gtf2gff.pl``); ``gff3=True`` writes GFF3 with ``ID``/``Parent`` attributes."""
    target = _out(output_path)
    arguments = [f"--out={target}"]
    for flag_name, value, label in (
        ("--gff3", gff3, "gff3"),
        ("--printExon", print_exons, "print_exons"),
        ("--printUTR", print_utrs, "print_utrs"),
        ("--printIntron", print_introns, "print_introns"),
        ("--includeStopInCDS", include_stop_in_cds, "include_stop_in_cds"),
    ):
        arguments += _a.switch(flag_name, _a.flag(value, label))
    return _a.simple(
        "gtf_to_gff",
        "gtf2gff.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        stdin=_in(gtf_path, "gtf_path"),
        to_stdout=False,
    )


def augustus_gth_to_gtf(
    gth_gff3_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert GenomeThreader GFF3 alignments to GTF (``gth2gtf.pl``)."""
    target = _out(output_path)
    return _a.simple(
        "gth_to_gtf",
        "gth2gtf.pl",
        [str(_in(gth_gff3_path, "gth_gff3_path")), str(target)],
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_mask_n_regions(
    genome_path: str, output_path: str, *, augustus_bin_dir: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Write a GFF-like list of the N (gap) regions of a genome FASTA (``maskNregions.pl``)."""
    return _a.simple(
        "mask_n_regions",
        "maskNregions.pl",
        [str(_in(genome_path, "genome_path"))],
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_random_split(
    genbank_path: str,
    test_size: int,
    train_path: str,
    test_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly split a GenBank gene set into a test set of ``test_size`` loci and a training set (``randomSplit.pl``)."""
    source = _in(genbank_path, "genbank_path")
    train = _out(train_path, "train_path")
    test = _out(test_path, "test_path")
    size = _cli.integer(test_size, "test_size", minimum=1)
    with tempfile.TemporaryDirectory(prefix="dnakit-augustus-split-") as scratch:
        copy = Path(scratch) / "genes.gb"
        shutil.copyfile(source, copy)
        completed = _a.run(
            "randomSplit.pl",
            [str(copy), str(size)],
            directory=_a.bin_dir(augustus_bin_dir),
            cwd=Path(scratch),
            capture=True,
            expected=[Path(str(copy) + ".train"), Path(str(copy) + ".test")],
        )
        _a.captured_text(completed)
        shutil.move(str(copy) + ".train", train)
        shutil.move(str(copy) + ".test", test)

    def loci(path: Path) -> int:
        return sum(
            1
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.startswith("LOCUS")
        )

    return _a.result(
        "random_split",
        completed,
        outputs={"train": train, "test": test},
        train_loci=loci(train),
        test_loci=loci(test),
    )


def augustus_summarize_base_content(
    fasta_path: str, output_path: str, *, augustus_bin_dir: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Summarise the A, C, G, T, N and other base counts of each sequence (``summarizeACGTcontent.pl``)."""
    return _a.simple(
        "summarize_base_content",
        "summarizeACGTcontent.pl",
        [str(_in(fasta_path, "fasta_path"))],
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


__all__ = [
    "augustus_add_names_to_gff3",
    "augustus_cegma_to_gff",
    "augustus_check_gff_names_in_fasta",
    "augustus_cluster_gene_dependencies",
    "augustus_compare_masking",
    "augustus_compute_flanking_region",
    "augustus_convert_genepred",
    "augustus_extract_annotation_fasta",
    "augustus_extract_cds_sequences",
    "augustus_extract_mrna_sequences",
    "augustus_filter_genbank_genes",
    "augustus_fix_joingenes_gtf",
    "augustus_gbrowse_old_to_gff3",
    "augustus_genbank_loci",
    "augustus_genbank_to_gff",
    "augustus_gff_to_genbank",
    "augustus_gtf_to_gff",
    "augustus_gth_to_gtf",
    "augustus_mask_n_regions",
    "augustus_merge_masking",
    "augustus_partition_gene_models",
    "augustus_random_split",
    "augustus_scipio_gff_to_gff",
    "augustus_scipio_yaml_to_gff",
    "augustus_summarize_base_content",
]
_register_public_verbose(globals())
