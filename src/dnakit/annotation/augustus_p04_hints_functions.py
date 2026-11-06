"""AUGUSTUS hint generation and alignment/coverage converters (integration pack P04, part 2).

Each public function runs one AUGUSTUS 3.5.0 auxiliary program or helper script (``bam2hints``, ``filterBam``,
``bam2wig``/``bamToWig.py``, ``aln2wig``, ``compileSpliceCands``, ``blat2hints.pl``, ``wig2hints.pl``, ``filterPSL.pl``
and others) with an explicit argument list and no shell. The user installs AUGUSTUS (``conda install -c bioconda
augustus``) and puts its ``bin`` folder on PATH or passes ``augustus_bin_dir``. Outputs are new explicit files;
nothing is overwritten. Many of these scripts need their input sorted (by query name, by target and start); the
facades do not sort. Results are JSON-serialisable dictionaries with status, backend, version, command, outputs,
small parsed counts and a bounded log tail. Defaults are those of the underlying programs.
"""

from __future__ import annotations

from collections.abc import Sequence
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


def _opts(**values: object) -> list[str]:
    """``--name=value`` for every value that is not None (False/True switches are handled by callers)."""
    return [f"--{key}={value}" for key, value in values.items() if value is not None]


def _i(value: object, name: str, minimum: int = 0) -> int | None:
    return _a.optional_int(value, name, minimum=minimum)


def _token(value: object, name: str) -> str | None:
    return None if value is None else _cli.text(value, name)


def augustus_bam_to_hints(
    bam_path: str,
    output_path: str,
    *,
    priority: int | None = None,
    max_gap_length: int | None = None,
    min_intron_length: int | None = None,
    max_intron_length: int | None = None,
    min_end_block_length: int | None = None,
    max_query_gap_length: int | None = None,
    exon_hints: bool = False,
    ep_cutoff: int | None = None,
    source: str | None = None,
    no_multiplicity: bool = False,
    remove_redundant: bool = False,
    max_coverage: int | None = None,
    splice_site_hints: bool = False,
    truncated_splice_site_hints: bool = False,
    score: float | None = None,
    max_gene_length: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Derive intron (and optionally exon) hints from spliced RNA-Seq alignments in BAM (``bam2hints``).

    The BAM must be sorted by target name and start. ``exon_hints`` adds exonpart, exon and splice-site hints;
    ``no_multiplicity`` keeps identical intron hints separate instead of summarising them with ``mult=n``.
    Unset options keep the program default (priority 4, max gap 14, min intron 32, max intron 350000, source E).
    """
    bam = _in(bam_path, "bam_path")
    target = _out(output_path)
    arguments = [f"--in={bam}", f"--out={target}"]
    arguments += _opts(
        priority=_i(priority, "priority"),
        maxgaplen=_i(max_gap_length, "max_gap_length"),
        minintronlen=_i(min_intron_length, "min_intron_length"),
        maxintronlen=_i(max_intron_length, "max_intron_length"),
        MinEndBlockLen=_i(min_end_block_length, "min_end_block_length"),
        maxQgaplen=_i(max_query_gap_length, "max_query_gap_length"),
        ep_cutoff=_i(ep_cutoff, "ep_cutoff"),
        source=_token(source, "source"),
        maxcoverage=_i(max_coverage, "max_coverage"),
        score=_a.optional_number(score, "score"),
        maxgenelen=_i(max_gene_length, "max_gene_length"),
    )
    for flag_name, value, label in (
        ("--exonhints", exon_hints, "exon_hints"),
        ("--nomult", no_multiplicity, "no_multiplicity"),
        ("--remove_redundant", remove_redundant, "remove_redundant"),
        ("--ssOn", splice_site_hints, "splice_site_hints"),
        ("--trunkSS", truncated_splice_site_hints, "truncated_splice_site_hints"),
    ):
        arguments += _a.switch(flag_name, _a.flag(value, label))
    return _a.simple(
        "bam_to_hints",
        "bam2hints",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_psl_to_hints(
    psl_path: str,
    output_path: str,
    *,
    priority: int | None = None,
    max_gap_length: int | None = None,
    min_intron_length: int | None = None,
    max_intron_length: int | None = None,
    max_query_gap_length: int | None = None,
    ep_cutoff: int | None = None,
    source: str | None = None,
    introns_only: bool = False,
    no_multiplicity: bool = False,
    remove_redundant: bool = False,
    max_coverage: int | None = None,
    splice_site_hints: bool = False,
    column_offset: int | None = None,
    score: float | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert BLAT or GMAP PSL cDNA alignments to AUGUSTUS hints (``blat2hints.pl``).

    The PSL must be sorted by target name (and start). ``column_offset=1`` reads the PSL layout of the UCSC database.
    Unset options keep the script defaults (priority 4, max gap 14, min intron 41, max intron 350000, source E).
    """
    psl = _in(psl_path, "psl_path")
    target = _out(output_path)
    arguments = [f"--in={psl}", f"--out={target}"]
    arguments += _opts(
        priority=_i(priority, "priority"),
        maxgaplen=_i(max_gap_length, "max_gap_length"),
        minintronlen=_i(min_intron_length, "min_intron_length"),
        maxintronlen=_i(max_intron_length, "max_intron_length"),
        maxQgaplen=_i(max_query_gap_length, "max_query_gap_length"),
        ep_cutoff=_i(ep_cutoff, "ep_cutoff"),
        source=_token(source, "source"),
        maxcoverage=_i(max_coverage, "max_coverage"),
        coloffset=_i(column_offset, "column_offset"),
        score=_a.optional_number(score, "score"),
    )
    for flag_name, value, label in (
        ("--intronsonly", introns_only, "introns_only"),
        ("--nomult", no_multiplicity, "no_multiplicity"),
        ("--remove_redundant", remove_redundant, "remove_redundant"),
        ("--ssOn", splice_site_hints, "splice_site_hints"),
    ):
        arguments += _a.switch(flag_name, _a.flag(value, label))
    return _a.simple(
        "psl_to_hints",
        "blat2hints.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_exonerate_to_hints(
    exonerate_path: str,
    output_path: str,
    *,
    priority: int | None = None,
    min_intron_length: int | None = None,
    max_intron_length: int | None = None,
    cds_part_cutoff: int | None = None,
    source: str | None = None,
    splice_site_hints: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert Exonerate protein2genome GFF alignments to AUGUSTUS hints (``exonerate2hints.pl``)."""
    source_file = _in(exonerate_path, "exonerate_path")
    target = _out(output_path)
    arguments = [f"--in={source_file}", f"--out={target}"]
    arguments += _opts(
        priority=_i(priority, "priority"),
        minintronlen=_i(min_intron_length, "min_intron_length"),
        maxintronlen=_i(max_intron_length, "max_intron_length"),
        CDSpart_cutoff=_i(cds_part_cutoff, "cds_part_cutoff"),
        source=_token(source, "source"),
    )
    arguments += _a.switch("--ssOn", _a.flag(splice_site_hints, "splice_site_hints"))
    return _a.simple(
        "exonerate_to_hints",
        "exonerate2hints.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def _stdin_stdout(
    operation: str,
    tool: str,
    input_path: str,
    output_path: str,
    arguments: list[str],
    directory: str | None,
    **extra: object,
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


def augustus_exoniphy_to_hints(
    table_path: str,
    output_path: str,
    *,
    priority: int | None = None,
    cutoff: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert exoniphy exon predictions (tab table: name, chrom, strand, cdsStart, cdsEnd) to hints (``exoniphyDb2hints.pl``)."""
    arguments = _opts(priority=_i(priority, "priority"), cutoff=_i(cutoff, "cutoff"))
    return _stdin_stdout(
        "exoniphy_to_hints",
        "exoniphyDb2hints.pl",
        table_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_polya_to_hints(
    input_path: str,
    output_path: str,
    *,
    input_format: str = "table",
    table_format: str | None = None,
    hint_radius: int | None = None,
    swap_strand: bool = False,
    radius: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert polyA/polyT positions to transcription-end (tts) hints.

    ``input_format="table"`` reads a polyA/polyT table (``polyA2hints.pl``; ``table_format`` ``new`` or ``newer``,
    ``hint_radius``, ``swap_strand`` for polyT files); ``input_format="pasa"`` reads PASA polyA output in FASTA format
    (``pasapolyA2hints.pl``, ``radius``).
    """
    source = _in(input_path, "input_path")
    target = _out(output_path)
    directory = _a.bin_dir(augustus_bin_dir)
    if _cli.choice(input_format, ("table", "pasa"), "input_format") == "table":
        arguments = [f"--in={source}", f"--out={target}"]
        if table_format is not None:
            arguments.append(
                f"--format={_cli.choice(table_format, ('new', 'newer'), 'table_format')}"
            )
        arguments += _opts(hintradius=_i(hint_radius, "hint_radius"))
        arguments += _a.switch("--swapstrand", _a.flag(swap_strand, "swap_strand"))
        return _a.simple(
            "polya_to_hints",
            "polyA2hints.pl",
            arguments,
            directory=directory,
            output=target,
            to_stdout=False,
        )
    arguments = [str(source)] + _opts(radius=_i(radius, "radius"))
    return _a.simple(
        "polya_to_hints", "pasapolyA2hints.pl", arguments, directory=directory, output=target
    )


def augustus_peptides_to_hints(
    psl_path: str,
    gff_path: str,
    fasta_path: str,
    output_path: str,
    *,
    source: str = "P",
    priority: int = 4,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert peptide-to-genome PSL alignments (with peptide GFF and FASTA) to hints (``peptides2hints.pl``)."""
    arguments = [
        str(_in(psl_path, "psl_path")),
        str(_in(gff_path, "gff_path")),
        str(_in(fasta_path, "fasta_path")),
        _cli.text(source, "source"),
        str(_cli.integer(priority, "priority", minimum=0)),
    ]
    return _a.simple(
        "peptides_to_hints",
        "peptides2hints.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_phastcons_to_hints(
    table_path: str,
    output_path: str,
    *,
    priority: int | None = None,
    source: str | None = None,
    min_score: int | None = None,
    cut_ends: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a phastCons conserved-element table to hints (``phastconsDB2hints.pl``)."""
    arguments = _opts(
        priority=_i(priority, "priority"),
        source=_token(source, "source"),
        minscore=_i(min_score, "min_score"),
        cutends=_i(cut_ends, "cut_ends"),
    )
    return _stdin_stdout(
        "phastcons_to_hints",
        "phastconsDB2hints.pl",
        table_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_retrogenes_to_hints(
    table_path: str,
    output_path: str,
    *,
    priority: int | None = None,
    source: str | None = None,
    min_score: int | None = None,
    max_small_gap_length: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a UCSC retroGenes table to hints (``retroDB2hints.pl``)."""
    arguments = _opts(
        priority=_i(priority, "priority"),
        source=_token(source, "source"),
        minscore=_i(min_score, "min_score"),
        maxSmallGapLen=_i(max_small_gap_length, "max_small_gap_length"),
    )
    return _stdin_stdout(
        "retrogenes_to_hints",
        "retroDB2hints.pl",
        table_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_transmap_to_hints(
    genepred_path: str,
    output_path: str,
    *,
    keep_ids_path: str | None = None,
    priority: int | None = None,
    ep_cutoff: int | None = None,
    ep_margin: int | None = None,
    ip_cutoff: int | None = None,
    utr_end_cutoff: int | None = None,
    min_intron_length: int | None = None,
    min_intron_length_utr: int | None = None,
    start_stop_radius: int | None = None,
    tss_tts_radius: int | None = None,
    source: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert transMap alignments (genePred-like table) to hints (``transMap2hints.pl``); ``keep_ids_path`` restricts the genes."""
    source_file = _in(genepred_path, "genepred_path")
    target = _out(output_path)
    arguments = [f"--in={source_file}", f"--out={target}"]
    if keep_ids_path is not None:
        arguments.append(f"--keepids={_in(keep_ids_path, 'keep_ids_path')}")
    arguments += _opts(
        priority=_i(priority, "priority"),
        ep_cutoff=_i(ep_cutoff, "ep_cutoff"),
        ep_margin=_i(ep_margin, "ep_margin"),
        ip_cutoff=_i(ip_cutoff, "ip_cutoff"),
        utrend_cutoff=_i(utr_end_cutoff, "utr_end_cutoff"),
        min_intron_len=_i(min_intron_length, "min_intron_length"),
        min_intron_len_utr=_i(min_intron_length_utr, "min_intron_length_utr"),
        start_stop_radius=_i(start_stop_radius, "start_stop_radius"),
        tss_tts_radius=_i(tss_tts_radius, "tss_tts_radius"),
        source=_token(source, "source"),
    )
    return _a.simple(
        "transmap_to_hints",
        "transMap2hints.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_wiggle_to_hints(
    wiggle_path: str,
    output_path: str,
    *,
    width: int | None = None,
    radius: float | None = None,
    margin: int | None = None,
    min_score: float | None = None,
    min_threshold: float | None = None,
    hint_type: str | None = None,
    prune: float | None = None,
    priority: int | None = None,
    strand: str | None = None,
    source: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a wiggle coverage (or conservation) track to exonpart/CDSpart hints (``wig2hints.pl``)."""
    arguments = _opts(
        width=_i(width, "width"),
        radius=_a.optional_number(radius, "radius"),
        margin=_i(margin, "margin"),
        minscore=_a.optional_number(min_score, "min_score"),
        minthresh=_a.optional_number(min_threshold, "min_threshold"),
        type=_token(hint_type, "hint_type"),
        prune=_a.optional_number(prune, "prune"),
        pri=_i(priority, "priority"),
        strand=_token(strand, "strand"),
        src=_token(source, "source"),
    )
    return _stdin_stdout(
        "wiggle_to_hints", "wig2hints.pl", wiggle_path, output_path, arguments, augustus_bin_dir
    )


def augustus_join_multiplicity_hints(
    hints_path: str, output_path: str, *, augustus_bin_dir: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Collapse identical hints of one source into one hint with ``mult=n`` (``join_mult_hints.pl``); the hints must be sorted."""
    return _stdin_stdout(
        "join_multiplicity_hints",
        "join_mult_hints.pl",
        hints_path,
        output_path,
        [],
        augustus_bin_dir,
    )


def augustus_remove_redundant_hints(
    hints_path: str,
    reference_hints_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove hints from ``hints_path`` that are redundant with ``reference_hints_path`` (``rmRedundantHints.pl``)."""
    arguments = [
        str(_in(hints_path, "hints_path")),
        str(_in(reference_hints_path, "reference_hints_path")),
    ]
    return _a.simple(
        "remove_redundant_hints",
        "rmRedundantHints.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_filter_splice_hints(
    genome_path: str,
    hints_path: str,
    output_path: str,
    *,
    splice_pattern: str = "GTAG",
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep intron hints whose splice sites match ``splice_pattern`` (usually GTAG, sometimes GCAG) in the genome (``filterSpliceHints.pl``)."""
    pattern = _cli.text(splice_pattern, "splice_pattern")
    if len(pattern) != 4 or not pattern.isalpha():
        raise ConfigurationError("splice_pattern must be four letters such as GTAG.")
    arguments = [str(_in(genome_path, "genome_path")), str(_in(hints_path, "hints_path")), pattern]
    return _a.simple(
        "filter_splice_hints",
        "filterSpliceHints.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_compile_splice_candidates(
    candidates_path: str,
    genome_path: str,
    output_path: str,
    *,
    compare_count: int | None = None,
    coverage_tolerance: float | None = None,
    margin: int | None = None,
    max_intron_length: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compile intron candidates from exon-candidate and coverage hints (``compileSpliceCands``); GFF goes to ``output_path``."""
    arguments = [
        "-f",
        str(_in(candidates_path, "candidates_path")),
        "-g",
        str(_in(genome_path, "genome_path")),
    ]
    arguments += _a.valued("-c", _i(compare_count, "compare_count", 1))
    arguments += _a.valued(
        "-t", _a.optional_number(coverage_tolerance, "coverage_tolerance", minimum=0, maximum=1)
    )
    arguments += _a.valued("-m", _i(margin, "margin", 1))
    arguments += _a.valued("-i", _i(max_intron_length, "max_intron_length", 1))
    return _a.simple(
        "compile_splice_candidates",
        "compileSpliceCands",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_choose_wiggle_intervals(
    wiggle_path: str,
    output_path: str,
    *,
    min_coverage: int | None = None,
    min_relative_coverage: float | None = None,
    max_gap: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Choose, per target sequence, an interval with good wiggle coverage and write it as BED (``wigchoose.pl``)."""
    arguments = _opts(
        mincov=_i(min_coverage, "min_coverage"),
        minrelcov=_a.optional_number(
            min_relative_coverage, "min_relative_coverage", minimum=0, maximum=1
        ),
        maxgap=_i(max_gap, "max_gap"),
    )
    return _stdin_stdout(
        "choose_wiggle_intervals",
        "wigchoose.pl",
        wiggle_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_bam_to_wiggle(
    bam_path: str,
    output_path: str,
    *,
    method: str = "bam2wig",
    genome_path: str | None = None,
    track_name: str | None = None,
    region: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert BAM alignments to wiggle coverage for exonpart hints.

    ``method="bam2wig"`` runs the AUGUSTUS program (sorted BAM; ``region`` such as ``chr1:10-250`` needs a BAM index;
    ``track_name``). ``method="bamToWig.py"`` runs the Python alternative, which also needs ``genome_path`` (FASTA) and
    samtools plus the UCSC tools ``twoBitInfo`` and ``faToTwoBit`` on PATH.
    """
    bam = _in(bam_path, "bam_path")
    target = _out(output_path)
    directory = _a.bin_dir(augustus_bin_dir)
    if _cli.choice(method, ("bam2wig", "bamToWig.py"), "method") == "bam2wig":
        arguments = []
        if region is not None:
            arguments += ["-r", _cli.text(region, "region")]
        arguments += _a.valued("-t", _token(track_name, "track_name"))
        arguments.append(str(bam))
        return _a.simple("bam_to_wiggle", "bam2wig", arguments, directory=directory, output=target)
    if genome_path is None:
        raise ConfigurationError("method='bamToWig.py' needs genome_path.")
    arguments = ["-b", str(bam), "-g", str(_in(genome_path, "genome_path")), "-o", str(target)]
    return _a.simple(
        "bam_to_wiggle",
        "bamToWig.py",
        arguments,
        directory=directory,
        output=target,
        to_stdout=False,
    )


def augustus_alignments_to_wiggle(
    alignment_path: str,
    output_path: str,
    *,
    span: bool = False,
    track_name: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert PSL or SHRIMP alignments to wiggle coverage (``aln2wig -f``); ``span`` uses span notation."""
    arguments = ["-f", str(_in(alignment_path, "alignment_path"))]
    arguments += _a.switch("-s", _a.flag(span, "span"))
    arguments += _a.valued("-t", _token(track_name, "track_name"))
    return _a.simple(
        "alignments_to_wiggle",
        "aln2wig",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_bedgraph_to_wiggle(
    bedgraph_path: str,
    output_path: str,
    *,
    span: bool = False,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert bedGraph coverage to variableStep wiggle (``bedgraph2wig.pl``); ``span`` uses span notation."""
    target = _out(output_path)
    arguments = [f"--bedgraphfile={_in(bedgraph_path, 'bedgraph_path')}", f"--outputfile={target}"]
    arguments += _a.switch("--span", _a.flag(span, "span"))
    return _a.simple(
        "bedgraph_to_wiggle",
        "bedgraph2wig.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_maf_to_conservation_wiggle(
    maf_path: str,
    output_path: str,
    *,
    window_size: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a MAF alignment to a conservation wiggle track (``maf2conswig.pl``)."""
    arguments = _opts(winsize=_i(window_size, "window_size", 1))
    return _stdin_stdout(
        "maf_to_conservation_wiggle",
        "maf2conswig.pl",
        maf_path,
        output_path,
        arguments,
        augustus_bin_dir,
    )


def augustus_filter_bam(
    bam_path: str,
    output_path: str,
    *,
    unique: bool = False,
    unique_threshold: float | None = None,
    best: bool = False,
    min_identity: int | None = None,
    min_coverage: int | None = None,
    no_introns: bool = False,
    insert_limit: int | None = None,
    paired: bool = False,
    max_intron_length: int | None = None,
    pair_bed_path: str | None = None,
    common_gene_path: str | None = None,
    threads: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter BAM alignments by identity, coverage, uniqueness and pairing (``filterBam``).

    The BAM must be sorted by query name. Program defaults: min identity 92, min coverage 80, unique threshold
    0.96. ``pair_bed_path`` and ``common_gene_path`` are new optional output files. Returns the record count of the
    filtered BAM when ``samtools`` is available.
    """
    bam = _in(bam_path, "bam_path")
    target = _out(output_path)
    arguments = ["--in", str(bam), "--out", str(target)]
    outputs_extra: dict[str, Path] = {}
    arguments += _a.switch("--uniq", _a.flag(unique, "unique")) + _a.switch(
        "--best", _a.flag(best, "best")
    )
    arguments += _a.valued(
        "--uniqThresh",
        _a.optional_number(unique_threshold, "unique_threshold", minimum=0, maximum=1),
    )
    arguments += _a.valued("--minId", _i(min_identity, "min_identity")) + _a.valued(
        "--minCover", _i(min_coverage, "min_coverage")
    )
    arguments += _a.switch("--noIntrons", _a.flag(no_introns, "no_introns")) + _a.valued(
        "--insertLimit", _i(insert_limit, "insert_limit")
    )
    arguments += _a.switch("--paired", _a.flag(paired, "paired")) + _a.valued(
        "--maxIntronLen", _i(max_intron_length, "max_intron_length")
    )
    if pair_bed_path is not None:
        outputs_extra["pair_bed"] = _out(pair_bed_path, "pair_bed_path")
        arguments += ["--pairBedFile", str(outputs_extra["pair_bed"])]
    if common_gene_path is not None:
        outputs_extra["common_genes"] = _out(common_gene_path, "common_gene_path")
        arguments += ["--commonGeneFile", str(outputs_extra["common_genes"])]
    arguments += _a.valued("--threads", _i(threads, "threads", 1))
    completed = _a.run(
        "filterBam",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=target.parent,
        capture=True,
        expected=[target],
    )
    _a.captured_text(completed)
    return _a.result(
        "filter_bam",
        completed,
        outputs={"output": target, **outputs_extra},
        bytes=target.stat().st_size,
    )


def _filter_stdin(
    operation: str,
    tool: str,
    input_path: str,
    output_path: str,
    arguments: list[str],
    directory: str | None,
    extras: dict[str, str | None],
) -> dict[str, Any]:
    target = _out(output_path)
    outputs: dict[str, Path] = {}
    for key, value in extras.items():
        if value is not None:
            outputs[key] = _out(value, f"{key}_path")
    for key, path in outputs.items():
        flag_name = {"pair_bed": "--pairbed", "common_genes": "--commongenefile"}[key]
        arguments.append(f"{flag_name}={path}")
    completed = _a.run(
        tool,
        arguments,
        directory=_a.bin_dir(directory),
        cwd=target.parent,
        stdin=_in(input_path, "input_path"),
        stdout=target,
    )
    return _a.result(
        operation, completed, outputs={"output": target, **outputs}, **_a.line_summary(target)
    )


def augustus_filter_psl(
    psl_path: str,
    output_path: str,
    *,
    min_identity: int | None = None,
    min_coverage: int | None = None,
    unique: bool = False,
    unique_threshold: float | None = None,
    best: bool = False,
    no_introns: bool = False,
    paired: bool = False,
    max_intron_length: int | None = None,
    pair_bed_path: str | None = None,
    common_gene_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter a PSL alignment file (BLAT, GMAP) by identity and coverage (``filterPSL.pl``); the PSL must be sorted by query name."""
    arguments = _opts(
        minId=_i(min_identity, "min_identity"),
        minCover=_i(min_coverage, "min_coverage"),
        uniqthresh=_a.optional_number(unique_threshold, "unique_threshold", minimum=0, maximum=1),
        maxintronlen=_i(max_intron_length, "max_intron_length"),
    )
    for flag_name, value, label in (
        ("--uniq", unique, "unique"),
        ("--best", best, "best"),
        ("--nointrons", no_introns, "no_introns"),
        ("--paired", paired, "paired"),
    ):
        arguments += _a.switch(flag_name, _a.flag(value, label))
    return _filter_stdin(
        "filter_psl",
        "filterPSL.pl",
        psl_path,
        output_path,
        arguments,
        augustus_bin_dir,
        {"pair_bed": pair_bed_path, "common_genes": common_gene_path},
    )


def augustus_filter_shrimp(
    shrimp_path: str,
    output_path: str,
    *,
    min_score: int | None = None,
    unique: bool = False,
    unique_threshold: float | None = None,
    best: bool = False,
    common_gene_path: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter a SHRIMP alignment file by score (``filterShrimp.pl``); the file must be sorted by query name."""
    arguments = _opts(
        minScore=_i(min_score, "min_score"),
        uniqthresh=_a.optional_number(unique_threshold, "unique_threshold", minimum=0, maximum=1),
    )
    arguments += _a.switch("--uniq", _a.flag(unique, "unique")) + _a.switch(
        "--best", _a.flag(best, "best")
    )
    return _filter_stdin(
        "filter_shrimp",
        "filterShrimp.pl",
        shrimp_path,
        output_path,
        arguments,
        augustus_bin_dir,
        {"common_genes": common_gene_path},
    )


def augustus_filter_maf(
    maf_path: str,
    output_path: str,
    *,
    species: Sequence[str] | None = None,
    min_sequences: int | None = None,
    interval: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep MAF blocks that intersect ``interval`` (``seqID:start-end``) or restrict them to ``species`` (``filterMaf.pl``)."""
    arguments: list[str] = []
    if species:
        for name in species:
            _a.name_token(name, "species")
        arguments.append("--species=" + ",".join(species))
    if min_sequences is not None:
        arguments += ["--min-seq", str(_i(min_sequences, "min_sequences", 2))]
    if interval is not None:
        text = _cli.text(interval, "interval")
        if ":" not in text or "-" not in text.split(":", 1)[1]:
            raise ConfigurationError("interval must look like seqID:start-end.")
        arguments.append(f"--interval={text}")
    return _stdin_stdout(
        "filter_maf", "filterMaf.pl", maf_path, output_path, arguments, augustus_bin_dir
    )


def augustus_weed_maf(
    maf_path: str, output_path: str, *, augustus_bin_dir: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Remove columns that contain only gaps from the blocks of a MAF alignment (``weedMaf.pl``)."""
    return _stdin_stdout("weed_maf", "weedMaf.pl", maf_path, output_path, [], augustus_bin_dir)


def augustus_split_hal_to_maf(
    hal_path: str,
    reference_genome: str,
    output_dir: str,
    *,
    reference_sequence: str | None = None,
    chunk_size: int | None = None,
    overlap: int | None = None,
    keep_dupes: bool = False,
    keep_ancestors: bool = False,
    no_split_list_path: str | None = None,
    cpus: int | None = None,
    hal_bin_dir: str | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a HAL genome alignment to MAF chunks of a given size with overlap (``hal2maf_split.pl``).

    Needs the HAL tools (``hal2maf``, ``halStats``) in ``hal_bin_dir`` or on PATH. ``no_split_list_path`` lists
    genic intervals (seqname, start, end) inside which the alignment may not be split. The MAF files are written to the
    new ``output_dir``.
    """
    hal = _in(hal_path, "hal_path")
    out = _cli.new_directory(output_dir, "output_dir")
    arguments = [
        "--halfile",
        str(hal),
        "--refGenome",
        _a.name_token(reference_genome, "reference_genome"),
        "--outdir",
        str(out),
    ]
    arguments += _a.valued(
        "--refSequence",
        None
        if reference_sequence is None
        else _a.name_token(reference_sequence, "reference_sequence"),
    )
    arguments += _a.valued("--chunksize", _i(chunk_size, "chunk_size", 1)) + _a.valued(
        "--overlap", _i(overlap, "overlap")
    )
    arguments += _a.switch("--keepDupes", _a.flag(keep_dupes, "keep_dupes")) + _a.switch(
        "--keepAncestors", _a.flag(keep_ancestors, "keep_ancestors")
    )
    if no_split_list_path is not None:
        arguments += ["--no_split_list", str(_in(no_split_list_path, "no_split_list_path"))]
    arguments += _a.valued("--cpus", _i(cpus, "cpus", 1))
    if hal_bin_dir is not None:
        folder = Path(hal_bin_dir).expanduser()
        if not folder.is_dir():
            raise ConfigurationError(f"hal_bin_dir is not a folder: {folder}")
        arguments += ["--hal_exec_dir", str(folder.resolve())]
    completed = _a.run(
        "hal2maf_split.pl", arguments, directory=_a.bin_dir(augustus_bin_dir), cwd=out
    )
    files = sorted(path.name for path in out.glob("*.maf"))
    return _a.result("split_hal_to_maf", completed, outputs={"output_dir": out}, maf_files=files)


def augustus_psl_map(
    psl_path: str,
    map_path: str,
    output_path: str,
    *,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Chain PSL alignments transitively through a second PSL map (``pslMap.pl --in --map``)."""
    target = _out(output_path)
    arguments = [
        f"--in={_in(psl_path, 'psl_path')}",
        f"--map={_in(map_path, 'map_path')}",
        f"--out={target}",
    ]
    return _a.simple(
        "psl_map",
        "pslMap.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=target,
        to_stdout=False,
    )


def augustus_map_junction_sam(
    sam_path: str,
    map_psl_path: str,
    output_path: str,
    *,
    flank: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map a SAM of reads aligned to an exon-exon junction file back to genome coordinates (``samMap.pl``); the default flank is 75 bp."""
    arguments = [str(_in(sam_path, "sam_path")), str(_in(map_psl_path, "map_psl_path"))]
    if flank is not None:
        arguments.append(str(_i(flank, "flank", 1)))
    return _a.simple(
        "map_junction_sam",
        "samMap.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        output=_out(output_path),
    )


def augustus_make_utr_training_set(
    codons_gff_path: str,
    genome_path: str,
    est_psl_path: str,
    output_prefix: str,
    *,
    only_best: bool = False,
    distance: int | None = None,
    augustus_bin_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a UTR training set (GenBank and GFF) from start/stop codons and filtered EST PSL alignments (``makeUtrTrainingSet.pl``).

    Writes ``<output_prefix>.gb`` and ``<output_prefix>.gff`` (both must not exist).
    """
    prefix = Path(output_prefix).expanduser().resolve()
    genbank = prefix.with_name(prefix.name + ".gb")
    gff = prefix.with_name(prefix.name + ".gff")
    for path in (genbank, gff):
        if path.exists():
            raise ConfigurationError(f"{path} already exists.")
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder for output_prefix does not exist: {prefix.parent}")
    arguments = _a.switch("--onlybest", _a.flag(only_best, "only_best")) + _opts(
        dist=_i(distance, "distance")
    )
    arguments += [
        str(_in(codons_gff_path, "codons_gff_path")),
        str(_in(genome_path, "genome_path")),
        str(_in(est_psl_path, "est_psl_path")),
        str(prefix),
    ]
    completed = _a.run(
        "makeUtrTrainingSet.pl",
        arguments,
        directory=_a.bin_dir(augustus_bin_dir),
        cwd=prefix.parent,
        capture=True,
        expected=[genbank, gff],
    )
    _a.captured_text(completed)
    return _a.result(
        "make_utr_training_set",
        completed,
        outputs={"genbank": genbank, "gff": gff},
        **_a.line_summary(gff),
    )


__all__ = [
    "augustus_alignments_to_wiggle",
    "augustus_bam_to_hints",
    "augustus_bam_to_wiggle",
    "augustus_bedgraph_to_wiggle",
    "augustus_choose_wiggle_intervals",
    "augustus_compile_splice_candidates",
    "augustus_exonerate_to_hints",
    "augustus_exoniphy_to_hints",
    "augustus_filter_bam",
    "augustus_filter_maf",
    "augustus_filter_psl",
    "augustus_filter_shrimp",
    "augustus_filter_splice_hints",
    "augustus_join_multiplicity_hints",
    "augustus_maf_to_conservation_wiggle",
    "augustus_make_utr_training_set",
    "augustus_map_junction_sam",
    "augustus_peptides_to_hints",
    "augustus_phastcons_to_hints",
    "augustus_polya_to_hints",
    "augustus_psl_map",
    "augustus_psl_to_hints",
    "augustus_remove_redundant_hints",
    "augustus_retrogenes_to_hints",
    "augustus_split_hal_to_maf",
    "augustus_transmap_to_hints",
    "augustus_weed_maf",
    "augustus_wiggle_to_hints",
]
_register_public_verbose(globals())
