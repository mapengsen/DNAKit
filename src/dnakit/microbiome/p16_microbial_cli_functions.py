"""Command-line facades for metagenome classification and microbial genome typing (pack P16).

Kaiju (kaijux, kaiju2krona, kaiju2table, kaiju-mergeOutputs, kaiju-addTaxonNames), Bracken (re-estimation report, k-mer
distribution build, table merge), geNomad, AMRFinderPlus (amrfinder, stxtyper, fasta_check, gff_check), ABRicate, mlst,
fastANI conserved-region mapping and plot, eggNOG-mapper, DIAMOND blastx and the two Prokka database converters.

Every function runs one program with a fixed argument whitelist, without a shell, through
``dnakit.bioinformatics._p16_cli``. A missing program raises ``BackendUnavailableError`` with an install hint. Results are
JSON-serialisable dictionaries (status, backend, version, command, outputs, log tail and small parsed tables). The
Bracken table merge and the fastANI plot are implemented in Python and need no program.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_KAIJU = cli.Tool(
    "kaiju", "kaiju", "Install Kaiju (conda install -c bioconda kaiju) and put it on PATH or pass the *_path argument.",
    "1.10.3", "https://doi.org/10.1038/ncomms11257", ("-h",), r"Kaiju (\d+(?:\.\d+)+)",
)
_BRACKEN = cli.Tool(
    "bracken", "bracken", "Install Bracken (conda install -c bioconda bracken kraken2) and put it on PATH or pass bracken_path.",
    "3.1", "https://doi.org/10.7717/peerj-cs.104", ("-v",), r"v?(\d+(?:\.\d+)+)",
)
_GENOMAD = cli.Tool(
    "genomad", "genomad", "Install geNomad (conda install -c bioconda genomad) and download its database "
    "(genomad download-database), then put genomad on PATH or pass genomad_path.",
    "1.11.0", "https://doi.org/10.1038/s41587-023-01953-y", ("--version",), None,
)
_AMRFINDER = cli.Tool(
    "amrfinder", "amrfinder", "Install AMRFinderPlus (conda install -c bioconda ncbi-amrfinderplus), run `amrfinder -u` "
    "once to fetch its database, then put it on PATH or pass the *_path argument.",
    "4.0", "https://doi.org/10.1038/s41598-021-91456-0", ("--version",), None,
)
_ABRICATE = cli.Tool(
    "abricate", "abricate", "Install ABRicate (conda install -c bioconda abricate) and put it on PATH or pass abricate_path.",
    "1.0.1", "https://github.com/tseemann/abricate", ("--version",), None,
)
_MLST = cli.Tool(
    "mlst", "mlst", "Install mlst (conda install -c bioconda mlst) and put it on PATH or pass mlst_path.",
    "2.23.0", "https://github.com/tseemann/mlst", ("--version",), None,
)
_FASTANI = cli.Tool(
    "fastani", "fastANI", "Install fastANI (conda install -c bioconda fastani) and put it on PATH or pass fastani_path.",
    "1.34", "https://doi.org/10.1038/s41467-018-07641-9", ("--version",), None,
)
_EMAPPER = cli.Tool(
    "eggnog-mapper", "emapper.py", "Install eggNOG-mapper (conda install -c bioconda eggnog-mapper), download its data "
    "(download_eggnog_data.py), then put emapper.py on PATH or pass emapper_path.",
    "2.1.12", "https://doi.org/10.1093/molbev/msab293", ("--version",), None,
)
_DIAMOND = cli.Tool(
    "diamond", "diamond", "Install DIAMOND (conda install -c bioconda diamond) and put it on PATH or pass diamond_path.",
    "2.2.8", "https://doi.org/10.1038/s41592-021-01101-x", ("version",), r"version (\d+(?:\.\d+)+)",
)
_PROKKA_DB = cli.Tool(
    "prokka", "prokka-genbank_to_fasta_db", "Install Prokka with BioPerl (conda install -c bioconda prokka) and put its "
    "bin folder on PATH or pass the *_path argument.",
    "1.14.6", "https://doi.org/10.1093/bioinformatics/btu153", ("--version",), None,
)
_RANKS = ("phylum", "class", "order", "family", "genus", "species")


def _output_file(value: str, name: str = "output_path") -> Path:
    return cli.new_file(value, name)


def _rank_list(value: Sequence[str] | None, name: str) -> str | None:
    if value is None:
        return None
    if isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of rank names.")
    return ",".join(cli.text(item, name) for item in value)


def _table(path: Path, names: Sequence[str] | None = None, header: bool = False, limit: int = 1000) -> list[dict[str, Any]]:
    return cli.read_table(path, names=names, header=header, limit=limit, text_columns=frozenset(n for n in (names or ()) if n != "score"))


# ---------------------------------------------------------------------------
# Kaiju
# ---------------------------------------------------------------------------


def search_kaiju_best_match(
    reads: str,
    database: str,
    output_path: str,
    *,
    reads_2: str | None = None,
    mode: str = "greedy",
    mismatches: int = 3,
    min_match_length: int = 11,
    min_score: int = 65,
    max_evalue: float = 0.01,
    seg_filter: bool = True,
    threads: int = 1,
    kaijux_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Best-matching database protein for DNA reads without taxonomy (``kaijux``).

    ``database`` is a Kaiju ``.fmi`` index of proteins; ``reads`` (and ``reads_2`` for pairs) are FASTA/FASTQ DNA reads.
    ``mode`` is ``greedy`` (uses ``mismatches``, ``min_score``, ``max_evalue``) or ``mem``. ``seg_filter`` is the SEG
    low-complexity filter (``-x``/``-X``). The output has classified (C/U), read name, score and the matched
    database identifiers per read; the first 1,000 rows are parsed into ``rows``.
    """
    source = cli.existing_file(reads, "reads")
    index = cli.existing_file(database, "database")
    target = _output_file(output_path)
    arguments = ["-f", str(index), "-i", str(source), "-o", str(target), "-z", str(cli.threads_value(threads))]
    if reads_2 is not None:
        arguments += ["-j", str(cli.existing_file(reads_2, "reads_2"))]
    arguments += ["-a", cli.choice(mode, ("greedy", "mem"), "mode")]
    arguments += ["-m", str(cli.integer(min_match_length, "min_match_length", minimum=1))]
    if mode == "greedy":
        arguments += ["-e", str(cli.integer(mismatches, "mismatches", minimum=0))]
        arguments += ["-s", str(cli.integer(min_score, "min_score", minimum=0))]
        arguments += ["-E", str(cli.number(max_evalue, "max_evalue", minimum=0.0))]
    arguments.append("-x" if cli.flag(seg_filter, "seg_filter") else "-X")
    record = cli.invoke(_KAIJU, arguments, cwd=target.parent, expected=[target], explicit=kaijux_path, executable_name="kaijux")
    rows = _table(target, ("classified", "read", "score", "matches"))
    return cli.finish(record, "search_kaiju_best_match", rows=rows, classified=sum(r["classified"] == "C" for r in rows))


def convert_kaiju_to_krona(
    kaiju_output: str,
    nodes: str,
    names: str,
    output_path: str,
    *,
    ranks: Sequence[str] | None = None,
    include_unclassified: bool = False,
    kaiju2krona_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert Kaiju output to a Krona text file (``kaiju2krona``).

    ``nodes``/``names`` are the taxonomy ``nodes.dmp``/``names.dmp``; ``ranks`` limits the taxon path to the listed ranks
    (``-l``); ``include_unclassified`` adds the count of unclassified reads (``-u``).
    """
    target = _output_file(output_path)
    arguments = ["-t", str(cli.existing_file(nodes, "nodes")), "-n", str(cli.existing_file(names, "names")),
                 "-i", str(cli.existing_file(kaiju_output, "kaiju_output")), "-o", str(target)]
    levels = _rank_list(ranks, "ranks")
    if levels:
        arguments += ["-l", levels]
    if cli.flag(include_unclassified, "include_unclassified"):
        arguments.append("-u")
    record = cli.invoke(_KAIJU, arguments, cwd=target.parent, expected=[target], explicit=kaiju2krona_path, executable_name="kaiju2krona")
    return cli.finish(record, "convert_kaiju_to_krona", line_count=cli.count_lines(target))


def summarize_kaiju_table(
    kaiju_outputs: Sequence[str],
    nodes: str,
    names: str,
    output_path: str,
    *,
    rank: str = "species",
    min_percent: float | None = None,
    min_reads: int | None = None,
    expand_viruses: bool = False,
    exclude_unclassified: bool = False,
    full_path: bool = False,
    path_ranks: Sequence[str] | None = None,
    kaiju2table_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-sample classification summary at a taxonomic rank (``kaiju2table``).

    ``rank`` is one of phylum, class, order, family, genus, species. Only one of ``min_percent`` (``-m``) and ``min_reads``
    (``-c``) may be given. ``expand_viruses`` (``-e``), ``exclude_unclassified`` (``-u``), ``full_path`` (``-p``) and
    ``path_ranks`` (``-l``) follow the program. One row per taxon is parsed into ``rows``.
    """
    sources = cli.existing_files(list(kaiju_outputs), "kaiju_outputs")
    target = _output_file(output_path)
    if min_percent is not None and min_reads is not None:
        raise ConfigurationError("Only one of min_percent and min_reads may be used.")
    arguments = ["-t", str(cli.existing_file(nodes, "nodes")), "-n", str(cli.existing_file(names, "names")),
                 "-r", cli.choice(rank, _RANKS, "rank"), "-o", str(target)]
    if min_percent is not None:
        arguments += ["-m", str(cli.number(min_percent, "min_percent", minimum=0.0, maximum=100.0))]
    if min_reads is not None:
        arguments += ["-c", str(cli.integer(min_reads, "min_reads", minimum=1))]
    for enabled, flag_name, switch in ((expand_viruses, "expand_viruses", "-e"), (exclude_unclassified, "exclude_unclassified", "-u"),
                                       (full_path, "full_path", "-p")):
        if cli.flag(enabled, flag_name):
            arguments.append(switch)
    levels = _rank_list(path_ranks, "path_ranks")
    if levels:
        arguments += ["-l", levels]
    arguments += [str(path) for path in sources]
    record = cli.invoke(_KAIJU, arguments, cwd=target.parent, expected=[target], explicit=kaiju2table_path, executable_name="kaiju2table")
    return cli.finish(record, "summarize_kaiju_table", rank=rank, rows=cli.read_table(target, limit=1000))


def merge_kaiju_outputs(
    first: str,
    second: str,
    output_path: str,
    *,
    conflict: str = "lca",
    nodes: str | None = None,
    use_score: bool = False,
    kaiju_merge_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the per-read classifications of two Kaiju/Kraken outputs (``kaiju-mergeOutputs``).

    Both inputs must be sorted by read name (column 2). ``conflict`` resolves disagreements: ``1`` or ``2`` (take that
    file), ``lca`` (default) or ``lowest``; ``lca`` and ``lowest`` need ``nodes`` (``nodes.dmp``). ``use_score`` prefers
    the better 4th-column score (``-s``).
    """
    target = _output_file(output_path)
    mode = cli.choice(conflict, ("1", "2", "lca", "lowest"), "conflict")
    arguments = ["-i", str(cli.existing_file(first, "first")), "-j", str(cli.existing_file(second, "second")),
                 "-o", str(target), "-c", mode]
    if mode in {"lca", "lowest"}:
        if nodes is None:
            raise ConfigurationError("conflict 'lca' and 'lowest' need nodes (nodes.dmp).")
        arguments += ["-t", str(cli.existing_file(nodes, "nodes"))]
    if cli.flag(use_score, "use_score"):
        arguments.append("-s")
    record = cli.invoke(_KAIJU, arguments, cwd=target.parent, expected=[target], explicit=kaiju_merge_path, executable_name="kaiju-mergeOutputs")
    return cli.finish(record, "merge_kaiju_outputs", conflict=mode, rows=_table(target, ("classified", "read", "taxon")))


def add_kaiju_taxon_names(
    kaiju_output: str,
    nodes: str,
    names: str,
    output_path: str,
    *,
    drop_unclassified: bool = False,
    full_path: bool = False,
    path_ranks: Sequence[str] | None = None,
    kaiju_names_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Append taxon names to Kaiju per-read output (``kaiju-addTaxonNames``).

    ``drop_unclassified`` (``-u``) removes unclassified reads, ``full_path`` (``-p``) prints the full lineage and
    ``path_ranks`` (``-r``) prints only the listed ranks.
    """
    target = _output_file(output_path)
    arguments = ["-t", str(cli.existing_file(nodes, "nodes")), "-n", str(cli.existing_file(names, "names")),
                 "-i", str(cli.existing_file(kaiju_output, "kaiju_output")), "-o", str(target)]
    if cli.flag(drop_unclassified, "drop_unclassified"):
        arguments.append("-u")
    if cli.flag(full_path, "full_path"):
        arguments.append("-p")
    levels = _rank_list(path_ranks, "path_ranks")
    if levels:
        arguments += ["-r", levels]
    record = cli.invoke(_KAIJU, arguments, cwd=target.parent, expected=[target], explicit=kaiju_names_path, executable_name="kaiju-addTaxonNames")
    return cli.finish(record, "add_kaiju_taxon_names", line_count=cli.count_lines(target))


# ---------------------------------------------------------------------------
# Bracken
# ---------------------------------------------------------------------------

_LEVELS = ("D", "P", "C", "O", "F", "G", "S")


def reestimate_bracken_report(
    kraken_report: str,
    database: str,
    output_table: str,
    output_report: str,
    *,
    read_length: int = 100,
    level: str = "S",
    threshold: int = 0,
    bracken_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bracken abundance table plus a Kraken-style report with the re-estimated read counts (``bracken -w``).

    ``kraken_report`` is a Kraken 2 report, ``database`` the Kraken database holding
    ``database<read_length>mers.kmer_distrib`` (see :func:`build_bracken_kmer_distribution`). ``level`` is D, P, C, O, F, G
    or S; ``threshold`` is the minimum read count before re-estimation (Bracken default 0).
    """
    source = cli.existing_file(kraken_report, "kraken_report")
    db = cli.existing_directory(database, "database")
    table = cli.new_file(output_table, "output_table")
    report = cli.new_file(output_report, "output_report")
    arguments = ["-d", str(db), "-i", str(source), "-o", str(table), "-w", str(report),
                 "-r", str(cli.integer(read_length, "read_length", minimum=1)), "-l", cli.choice(level, _LEVELS, "level"),
                 "-t", str(cli.integer(threshold, "threshold", minimum=0))]
    record = cli.invoke(_BRACKEN, arguments, cwd=table.parent, expected=[table, report], explicit=bracken_path)
    return cli.finish(record, "reestimate_bracken_report", rows=cli.read_table(table, limit=1000))


def build_bracken_kmer_distribution(
    database: str,
    *,
    read_length: int = 100,
    kmer_length: int = 35,
    threads: int = 1,
    kraken_type: str = "kraken2",
    bracken_build_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the read-length k-mer distribution of a Kraken database (``bracken-build``).

    Runs the two Bracken stages: simulate read-length k-mer assignments from the library (``kmer2read_distr`` ->
    ``database<L>mers.kraken``) and aggregate them (``generate_kmer_distribution.py`` -> ``database<L>mers.kmer_distrib``).
    Both files are written inside the existing Kraken ``database`` folder, as the program does; existing files are not
    overwritten (the call is refused). ``kraken_type`` is kraken2, kraken or krakenuniq.
    """
    db = cli.existing_directory(database, "database")
    length = cli.integer(read_length, "read_length", minimum=1)
    expected = [db / f"database{length}mers.kraken", db / f"database{length}mers.kmer_distrib"]
    for path in expected:
        if path.exists():
            raise ConfigurationError(f"Output already exists and will not be overwritten: {path}")
    arguments = ["-d", str(db), "-l", str(length), "-k", str(cli.integer(kmer_length, "kmer_length", minimum=1)),
                 "-t", str(cli.threads_value(threads)), "-y", cli.choice(kraken_type, ("kraken2", "kraken", "krakenuniq"), "kraken_type")]
    record = cli.invoke(_BRACKEN, arguments, cwd=db, expected=expected, explicit=bracken_build_path, executable_name="bracken-build")
    return cli.finish(record, "build_bracken_kmer_distribution", read_length=length)


def combine_bracken_tables(
    tables: Sequence[str],
    output_path: str,
    *,
    names: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge per-sample Bracken output tables into one multi-sample table (Bracken ``combine_bracken_outputs.py``).

    Implemented in Python with the script's semantics. Inputs are Bracken *output* tables (not reports) estimated at the
    same level. Columns: ``name``, ``taxonomy_id``, ``taxonomy_lvl`` and for each sample ``<sample>_num`` and
    ``<sample>_frac`` (fraction of the sample's re-estimated reads, five decimals). ``names`` give the sample names
    (default: the file base names).
    """
    sources = cli.existing_files(list(tables), "tables")
    target = cli.new_file(output_path, "output_path")
    samples = [cli.text(n, "names") for n in names] if names is not None else [p.name for p in sources]
    if len(samples) != len(sources) or len(set(samples)) != len(samples):
        raise ConfigurationError("names must be unique and match the number of tables.")
    counts: dict[str, dict[str, dict[str, int]]] = {}
    totals = {name: 0 for name in samples}
    level = ""
    for sample, path in zip(samples, sources, strict=True):
        with path.open("r", encoding="utf-8") as handle:
            next(handle, None)
            for line in handle:
                if not line.strip():
                    continue
                fields = line.strip().split("\t")
                if len(fields) != 7:
                    raise ConfigurationError(f"{path.name}: expected 7 tab-separated columns, got {len(fields)}.")
                name, taxid, taxlevel, _kraken, _added, estimated, _fraction = fields
                known = counts.setdefault(name, {})
                if known and taxid not in known:
                    raise ConfigurationError(f"Taxonomy IDs do not match for {name!r}.")
                if level and level != taxlevel:
                    raise ConfigurationError("Taxonomy level does not match between samples.")
                level = taxlevel
                totals[sample] += int(estimated)
                known.setdefault(taxid, {})[sample] = int(estimated)
    with target.open("x", encoding="utf-8") as out:
        out.write("name\ttaxonomy_id\ttaxonomy_lvl" + "".join(f"\t{s}_num\t{s}_frac" for s in samples) + "\n")
        for name, by_id in counts.items():
            taxid = next(iter(by_id))
            row = [name, taxid, level]
            for sample in samples:
                if sample in by_id[taxid]:
                    num = by_id[taxid][sample]
                    row += [str(num), f"{num / totals[sample]:0.5f}"]
                else:
                    row += ["0", "0.00000"]
            out.write("\t".join(row) + "\n")
    return {
        "status": "completed", "backend": "bracken-scripts", "operation": "combine_bracken_tables",
        "samples": samples, "taxa": len(counts), "outputs": [str(target)],
        "provenance": {"tool": "Bracken analysis_scripts/combine_bracken_outputs.py", "reimplemented_in": "dnakit"},
    }


# ---------------------------------------------------------------------------
# geNomad
# ---------------------------------------------------------------------------

_GENOMAD_STEPS = (
    "end-to-end", "annotate", "find-proviruses", "marker-classification", "nn-classification",
    "aggregated-classification", "score-calibration", "summary",
)
_GENOMAD_NEEDS_DB = {"end-to-end", "annotate", "find-proviruses", "marker-classification"}


def run_genomad(
    contigs: str,
    output_directory: str,
    *,
    step: str = "end-to-end",
    database: str | None = None,
    threads: int = 1,
    cleanup: bool = True,
    splits: int | None = None,
    min_score: float | None = None,
    enable_score_calibration: bool = False,
    genomad_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify virus and plasmid contigs with geNomad, whole or one stage at a time.

    ``step`` selects the sub-command: ``end-to-end`` (aggregated classification of the full pipeline), ``annotate``
    (gene prediction and marker annotation), ``find-proviruses``, ``marker-classification``, ``nn-classification``,
    ``aggregated-classification``, ``score-calibration`` (batch-corrected scores) and ``summary`` (filtering and
    reports). Stages after ``annotate`` read the intermediate files that earlier stages wrote into the same
    ``output_directory``; steps that need the geNomad database require ``database`` (the folder from
    ``genomad download-database``). ``output_directory`` is created if it does not exist and reused across stages.
    ``min_score`` is ``--min-score`` (summary/end-to-end), ``enable_score_calibration`` is
    ``--enable-score-calibration``, ``splits`` is ``--splits``.
    """
    fasta = cli.existing_file(contigs, "contigs")
    chosen = cli.choice(step, _GENOMAD_STEPS, "step")
    out_dir = Path(output_directory).expanduser().resolve()
    if out_dir.exists() and not out_dir.is_dir():
        raise ConfigurationError(f"output_directory is not a folder: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    arguments = [chosen]
    if cli.flag(cleanup, "cleanup") and chosen in {"end-to-end", "annotate"}:
        arguments.append("--cleanup")
    arguments += ["--threads", str(cli.threads_value(threads))]
    if splits is not None:
        arguments += ["--splits", str(cli.integer(splits, "splits", minimum=1))]
    if min_score is not None and chosen in {"end-to-end", "summary"}:
        arguments += ["--min-score", str(cli.number(min_score, "min_score", minimum=0.0, maximum=1.0))]
    if cli.flag(enable_score_calibration, "enable_score_calibration") and chosen in {"end-to-end", "aggregated-classification"}:
        arguments.append("--enable-score-calibration")
    arguments += [str(fasta), str(out_dir)]
    if chosen in _GENOMAD_NEEDS_DB:
        if database is None:
            raise ConfigurationError(f"step {chosen!r} needs database (the geNomad database folder).")
        arguments.append(str(cli.existing_directory(database, "database")))
    record = cli.invoke(_GENOMAD, arguments, cwd=out_dir, explicit=genomad_path)
    produced = sorted(str(p) for p in out_dir.rglob("*") if p.is_file())[:200]
    return cli.finish(record, "run_genomad", step=chosen, output_directory=str(out_dir), files=produced)


# ---------------------------------------------------------------------------
# AMRFinderPlus
# ---------------------------------------------------------------------------


def detect_amr_genes_nucleotide(
    contigs: str,
    output_path: str,
    *,
    organism: str | None = None,
    database: str | None = None,
    plus: bool = False,
    min_identity: float | None = None,
    min_coverage: float | None = None,
    translation_table: int = 11,
    name: str | None = None,
    threads: int = 1,
    amrfinder_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """AMR, stress and virulence genes and point mutations in nucleotide contigs (``amrfinder --nucleotide``).

    Without ``organism`` genes are reported (``amr_report``); with ``organism`` (see ``amrfinder --list_organisms``)
    organism-specific point mutations are also detected (``dna_mutation``). ``plus`` adds the plus genes; ``min_identity``
    and ``min_coverage`` are 0-1 (program defaults: curated threshold and 0.5). ``database`` is an AMRFinderPlus database
    folder (default: the installed one). The tab-separated report is written to ``output_path`` and parsed into ``rows``.
    """
    fasta = cli.existing_file(contigs, "contigs")
    target = cli.new_file(output_path, "output_path")
    arguments = ["--nucleotide", str(fasta), "--output", str(target), "--translation_table",
                 str(cli.integer(translation_table, "translation_table", minimum=1, maximum=33)),
                 "--threads", str(cli.threads_value(threads))]
    if organism is not None:
        arguments += ["--organism", cli.text(organism, "organism")]
    if database is not None:
        arguments += ["--database", str(cli.existing_directory(database, "database"))]
    if cli.flag(plus, "plus"):
        arguments.append("--plus")
    if min_identity is not None:
        arguments += ["--ident_min", str(cli.number(min_identity, "min_identity", minimum=-1.0, maximum=1.0))]
    if min_coverage is not None:
        arguments += ["--coverage_min", str(cli.number(min_coverage, "min_coverage", minimum=0.0, maximum=1.0))]
    if name is not None:
        arguments += ["--name", cli.text(name, "name")]
    record = cli.invoke(_AMRFINDER, arguments, cwd=target.parent, expected=[target], explicit=amrfinder_path)
    return cli.finish(record, "detect_amr_genes_nucleotide", organism=organism, rows=cli.read_table(target, limit=1000))


def type_shiga_toxin_operons(
    contigs: str,
    output_path: str,
    *,
    amrfinder_format: bool = False,
    name: str | None = None,
    threads: int = 1,
    stxtyper_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shiga-toxin (stx) operon type of nucleotide assemblies (``stxtyper``)."""
    fasta = cli.existing_file(contigs, "contigs")
    target = cli.new_file(output_path, "output_path")
    arguments = ["--nucleotide", str(fasta), "--output", str(target), "--threads", str(cli.threads_value(threads))]
    if cli.flag(amrfinder_format, "amrfinder_format"):
        arguments.append("--amrfinder")
    if name is not None:
        arguments += ["--name", cli.text(name, "name")]
    record = cli.invoke(_AMRFINDER, arguments, cwd=target.parent, expected=[target], explicit=stxtyper_path, executable_name="stxtyper")
    return cli.finish(record, "type_shiga_toxin_operons", rows=cli.read_table(target, limit=1000))


def check_nucleotide_fasta_amrfinder(
    fasta: str,
    *,
    fasta_check_path: str | None = None,
    work_directory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate a nucleotide FASTA and report the sequence count and lengths (AMRFinderPlus ``fasta_check``).

    ``fasta_check`` prints the number of sequences followed by the cumulative sequence length after each record; these
    numbers are parsed into ``sequence_count``, ``cumulative_lengths`` and the per-sequence ``lengths``. A malformed FASTA makes the program exit with an
    error, which is raised.
    """
    path = cli.existing_file(fasta, "fasta")
    cwd = path.parent if work_directory is None else cli.existing_directory(work_directory, "work_directory")
    record = cli.invoke(_AMRFINDER, [str(path)], cwd=cwd, explicit=fasta_check_path, executable_name="fasta_check")
    numbers = [int(token) for token in record["log_tail"].split() if token.isdigit()]
    cumulative = numbers[1:]
    lengths = [end - start for start, end in zip([0, *cumulative], cumulative, strict=False)]
    return cli.finish(record, "check_nucleotide_fasta_amrfinder", fasta=str(path),
                      sequence_count=numbers[0] if numbers else None, lengths=lengths, cumulative_lengths=cumulative)


def check_gff_annotation_amrfinder(
    gff: str,
    protein_fasta: str,
    nucleotide_fasta: str,
    *,
    annotation_format: str = "genbank",
    gff_check_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate a gene-annotation GFF against nucleotide and protein FASTA (AMRFinderPlus ``gff_check``).

    ``annotation_format`` is one of bakta, genbank, microscope, patric, pgap, prodigal, prokka, pseudomonasdb, rast,
    standard. Exit status 0 means the three files are consistent.
    """
    gff_path = cli.existing_file(gff, "gff")
    arguments = [str(gff_path), str(cli.existing_file(protein_fasta, "protein_fasta")),
                 str(cli.existing_file(nucleotide_fasta, "nucleotide_fasta")),
                 cli.choice(annotation_format, ("bakta", "genbank", "microscope", "patric", "pgap", "prodigal", "prokka",
                                                "pseudomonasdb", "rast", "standard"), "annotation_format")]
    record = cli.invoke(_AMRFINDER, arguments, cwd=gff_path.parent, explicit=gff_check_path, executable_name="gff_check")
    return cli.finish(record, "check_gff_annotation_amrfinder", gff=str(gff_path))


# ---------------------------------------------------------------------------
# ABRicate
# ---------------------------------------------------------------------------


def screen_contigs_abricate(
    contigs: str,
    output_path: str,
    *,
    database: str = "resfinder",
    min_identity: float = 75.0,
    min_coverage: float = 0.0,
    data_directory: str | None = None,
    abricate_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Screen contigs for AMR, virulence or plasmid genes by BLASTn (``abricate``).

    ``database`` is an installed ABRicate database (argannot, card, ncbi, ncbibetalactamase, plasmidfinder, resfinder,
    vfdb or a custom one); ``min_identity`` and ``min_coverage`` are percentages (ABRicate defaults 75 and 0).
    The tab-separated report is written to ``output_path`` and parsed into ``rows``.
    """
    fasta = cli.existing_file(contigs, "contigs")
    target = cli.new_file(output_path, "output_path")
    arguments = ["--db", cli.text(database, "database"), "--minid", str(cli.number(min_identity, "min_identity", minimum=0.0, maximum=100.0)),
                 "--mincov", str(cli.number(min_coverage, "min_coverage", minimum=0.0, maximum=100.0)), "--nopath"]
    if data_directory is not None:
        arguments += ["--datadir", str(cli.existing_directory(data_directory, "data_directory"))]
    arguments.append(str(fasta))
    record = cli.invoke(_ABRICATE, arguments, cwd=target.parent, stdout_to=target, explicit=abricate_path)
    return cli.finish(record, "screen_contigs_abricate", database=database, rows=cli.read_table(target, comment="#", limit=1000))


def summarize_abricate_reports(
    reports: Sequence[str],
    output_path: str,
    *,
    abricate_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene presence/absence matrix from several ABRicate reports (``abricate --summary``)."""
    sources = cli.existing_files(list(reports), "reports")
    target = cli.new_file(output_path, "output_path")
    record = cli.invoke(_ABRICATE, ["--summary", *map(str, sources)], cwd=target.parent, stdout_to=target, explicit=abricate_path)
    return cli.finish(record, "summarize_abricate_reports", rows=cli.read_table(target, limit=1000))


def download_abricate_database(
    database: str,
    *,
    database_directory: str | None = None,
    force: bool = False,
    abricate_get_db_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and format an ABRicate reference database (``abricate-get_db``; network operation).

    ``database`` is argannot, card, ncbi, ncbibetalactamase, plasmidfinder, resfinder or vfdb. The files are written under
    ``database_directory`` (default: the installation's ``db`` folder). Needs network access, BLAST+ and the program.
    """
    name = cli.choice(database, ("argannot", "card", "ncbi", "ncbibetalactamase", "plasmidfinder", "resfinder", "vfdb"), "database")
    arguments = ["--db", name]
    cwd = Path.cwd()
    if database_directory is not None:
        directory = Path(database_directory).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        arguments += ["--dbdir", str(directory)]
        cwd = directory
    if cli.flag(force, "force"):
        arguments.append("--force")
    record = cli.invoke(_ABRICATE, arguments, cwd=cwd, explicit=abricate_get_db_path, executable_name="abricate-get_db")
    return cli.finish(record, "download_abricate_database", database=name)


# ---------------------------------------------------------------------------
# mlst
# ---------------------------------------------------------------------------


def type_sequence_type_mlst(
    contigs: str,
    output_path: str,
    *,
    scheme: str | None = None,
    novel_alleles_path: str | None = None,
    min_identity: float = 95.0,
    min_coverage: float = 10.0,
    min_score: float = 50.0,
    threads: int = 1,
    mlst_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multi-locus sequence type of assembled contigs (``mlst``).

    The scheme is detected automatically unless ``scheme`` forces one. ``novel_alleles_path`` (a new FASTA file) receives
    novel alleles (``--novel``). Output columns: file, scheme, ST, then the allele calls (``rows`` holds the parsed
    first row as ``fields``). Scoring defaults follow mlst (95 / 10 / 50).
    """
    fasta = cli.existing_file(contigs, "contigs")
    target = cli.new_file(output_path, "output_path")
    arguments = ["--minid", str(cli.number(min_identity, "min_identity", minimum=0.0, maximum=100.0)),
                 "--mincov", str(cli.number(min_coverage, "min_coverage", minimum=0.0, maximum=100.0)),
                 "--minscore", str(cli.number(min_score, "min_score", minimum=0.0, maximum=100.0)),
                 "--threads", str(cli.threads_value(threads)), "--nopath"]
    if scheme is not None:
        arguments += ["--scheme", cli.text(scheme, "scheme")]
    novel = None
    if novel_alleles_path is not None:
        novel = cli.new_file(novel_alleles_path, "novel_alleles_path")
        arguments += ["--novel", str(novel)]
    arguments.append(str(fasta))
    record = cli.invoke(_MLST, arguments, cwd=target.parent, stdout_to=target, explicit=mlst_path)
    lines = [ln.split("\t") for ln in target.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if novel is not None and novel.exists():
        record["outputs"].append(str(novel))
    return cli.finish(record, "type_sequence_type_mlst", fields=lines[0] if lines else [])


# ---------------------------------------------------------------------------
# fastANI
# ---------------------------------------------------------------------------


def map_fastani_conserved_regions(
    query: str,
    reference: str,
    output_path: str,
    *,
    kmer: int = 16,
    fragment_length: int = 3000,
    min_fraction: float = 0.2,
    threads: int = 1,
    fastani_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """ANI of one genome pair plus the reciprocal conserved-region mapping records (``fastANI --visualize``).

    The summary line is returned (``ani``, ``matched_fragments``, ``total_fragments``); the per-fragment mapping file is
    ``<output_path>.visual`` and holds query/reference names, positions, identity and mapping coordinates.
    """
    q = cli.existing_file(query, "query")
    r = cli.existing_file(reference, "reference")
    target = cli.new_file(output_path, "output_path")
    visual = target.with_name(target.name + ".visual")
    if visual.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {visual}")
    arguments = ["-q", str(q), "-r", str(r), "-o", str(target), "--visualize", "-k", str(cli.integer(kmer, "kmer", minimum=1, maximum=16)),
                 "--fragLen", str(cli.integer(fragment_length, "fragment_length", minimum=1)),
                 "--minFraction", str(cli.number(min_fraction, "min_fraction", minimum=0.0, maximum=1.0)),
                 "-t", str(cli.threads_value(threads))]
    record = cli.invoke(_FASTANI, arguments, cwd=target.parent, expected=[target, visual], explicit=fastani_path)
    summary = target.read_text(encoding="utf-8").split()
    parsed: dict[str, Any] = {}
    if len(summary) >= 5:
        parsed = {"ani": float(summary[2]), "matched_fragments": int(summary[3]), "total_fragments": int(summary[4])}
    return cli.finish(record, "map_fastani_conserved_regions", visual_file=str(visual), **parsed)


def plot_fastani_conserved_regions(
    visual_file: str,
    image_path: str,
    *,
    min_identity: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dot plot of conserved regions between two genomes from a fastANI ``.visual`` file.

    The original R script (``scripts/visualize.R``) is not run; the same data (reference start/end against query
    start/end of every mapped fragment, coloured by identity) is drawn with matplotlib (``dnakit[bio-plots]``). Returns
    the PNG path and the fragment table (first 1,000 rows).
    """
    source = cli.existing_file(visual_file, "visual_file")
    image = cli.new_file(image_path, "image_path")
    plt = require("matplotlib.pyplot", "bio-plots")
    rows = cli.read_table(source, header=False, names=["query", "reference", "identity", "count_1", "count_2", "count_3", "query_start", "query_end", "ref_start", "ref_end", "extra_1", "extra_2"], text_columns=frozenset({"query", "reference", "count_1", "count_2", "count_3", "extra_1", "extra_2"}), limit=10**7)
    rows = [row for row in rows if isinstance(row["identity"], (int, float)) and row["identity"] >= min_identity]
    if not rows:
        raise ConfigurationError("The visual file has no mapped fragments at this identity.")
    figure, axis = plt.subplots(figsize=(6, 6))
    scatter = None
    for row in rows:
        scatter = axis.plot([row["ref_start"], row["ref_end"]], [row["query_start"], row["query_end"]], linewidth=1.5, color="#1f6f8b")
    axis.set_xlabel("reference position (bp)")
    axis.set_ylabel("query position (bp)")
    axis.set_title(f"{len(rows)} conserved fragments")
    figure.savefig(image, dpi=150, bbox_inches="tight")
    plt.close(figure)
    del scatter
    return {
        "status": "completed", "backend": "matplotlib", "operation": "plot_fastani_conserved_regions",
        "image": str(image), "fragments": len(rows), "rows": rows[:1000], "outputs": [str(image)],
        "provenance": {"tool": "fastANI scripts/visualize.R", "reimplemented_in": "dnakit"},
    }


# ---------------------------------------------------------------------------
# eggNOG-mapper and DIAMOND
# ---------------------------------------------------------------------------


def annotate_dna_eggnog_mapper(
    fasta: str,
    output_directory: str,
    *,
    prefix: str = "emapper",
    input_type: str = "CDS",
    gene_prediction: str | None = None,
    translate: bool = False,
    decorate_gff: str | None = None,
    data_directory: str | None = None,
    threads: int = 1,
    emapper_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Orthology annotation of DNA with eggNOG-mapper (``emapper.py``).

    ``input_type`` is ``CDS`` (coding sequences, optionally ``translate``), ``genome`` or ``metagenome``; for the last two
    genes are predicted first with ``gene_prediction`` = ``prodigal`` (ab initio) or ``search`` (homology-based blastx
    search), the program default when unset. ``decorate_gff`` is ``yes`` (create a GFF from the predictions) or the path of
    an existing GFF to decorate. ``data_directory`` is the eggNOG data folder (``--data_dir``). Needs the eggNOG
    database (about 50 GB with DIAMOND data); outputs are ``<prefix>.emapper.annotations`` (parsed into ``rows``).
    """
    source = cli.existing_file(fasta, "fasta")
    out_dir = cli.new_directory(output_directory, "output_directory")
    kind = cli.choice(input_type, ("CDS", "genome", "metagenome"), "input_type")
    arguments = ["-i", str(source), "--itype", kind, "--output_dir", str(out_dir), "-o", cli.text(prefix, "prefix"),
                 "--cpu", str(cli.threads_value(threads))]
    if gene_prediction is not None:
        if kind == "CDS":
            raise ConfigurationError("gene_prediction only applies to genome and metagenome input.")
        arguments += ["--genepred", cli.choice(gene_prediction, ("prodigal", "search"), "gene_prediction")]
    if cli.flag(translate, "translate"):
        arguments.append("--translate")
    if decorate_gff is not None:
        arguments += ["--decorate_gff", "yes" if decorate_gff == "yes" else str(cli.existing_file(decorate_gff, "decorate_gff"))]
    if data_directory is not None:
        arguments += ["--data_dir", str(cli.existing_directory(data_directory, "data_directory"))]
    annotations = out_dir / f"{prefix}.emapper.annotations"
    record = cli.invoke(_EMAPPER, arguments, cwd=out_dir, expected=[annotations], explicit=emapper_path)
    rows = cli.read_table(annotations, comment="#", limit=1000) if annotations.stat().st_size else []
    return cli.finish(record, "annotate_dna_eggnog_mapper", input_type=kind, gene_prediction=gene_prediction,
                      files=sorted(str(p) for p in out_dir.iterdir()), rows=rows)


def search_dna_translated_diamond(
    queries: str,
    database: str,
    output_path: str,
    *,
    output_format: str = "6",
    taxonomic_classification: bool = False,
    frameshift_penalty: int | None = None,
    long_reads: bool = False,
    exhaustive: bool = False,
    max_target_seqs: int | None = None,
    evalue: float = 0.001,
    sensitivity: str | None = None,
    threads: int = 1,
    diamond_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Translated search of DNA queries against a protein database (``diamond blastx``).

    ``database`` is a DIAMOND ``.dmnd`` file. ``taxonomic_classification=True`` writes per-read LCA taxonomy (output
    format 102; the database needs taxonomy data); otherwise ``output_format`` is a DIAMOND format number or name (``6``
    tabular default, ``0``, ``5``, ``100``...). ``frameshift_penalty`` (``--frameshift``) and ``long_reads``
    (``--long-reads``, exclusive with ``frameshift_penalty``) make the alignment frameshift-aware; ``exhaustive``
    (``--swipe``) searches all database sequences. ``sensitivity`` is one of fast, mid-sensitive, sensitive,
    more-sensitive, very-sensitive, ultra-sensitive. The first 1,000 rows of tabular output are parsed.
    """
    query = cli.existing_file(queries, "queries")
    db = cli.existing_file(database, "database")
    target = cli.new_file(output_path, "output_path")
    fmt = "102" if cli.flag(taxonomic_classification, "taxonomic_classification") else cli.text(str(output_format), "output_format")
    arguments = ["blastx", "--db", str(db), "--query", str(query), "--out", str(target), "--outfmt", fmt,
                 "--evalue", str(cli.number(evalue, "evalue", minimum=0.0)), "--threads", str(cli.threads_value(threads)), "--quiet"]
    if frameshift_penalty is not None and long_reads:
        raise ConfigurationError("frameshift_penalty and long_reads are mutually exclusive.")
    if frameshift_penalty is not None:
        arguments += ["--frameshift", str(cli.integer(frameshift_penalty, "frameshift_penalty", minimum=1))]
    if cli.flag(long_reads, "long_reads"):
        arguments.append("--long-reads")
    if cli.flag(exhaustive, "exhaustive"):
        arguments.append("--swipe")
    if max_target_seqs is not None:
        arguments += ["--max-target-seqs", str(cli.integer(max_target_seqs, "max_target_seqs", minimum=0))]
    if sensitivity is not None:
        arguments.append("--" + cli.choice(sensitivity, ("fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"), "sensitivity"))
    record = cli.invoke(_DIAMOND, arguments, cwd=target.parent, expected=[target], explicit=diamond_path)
    rows = cli.read_table(target, header=False, limit=1000) if fmt in {"6", "102"} and target.stat().st_size else []
    return cli.finish(record, "search_dna_translated_diamond", output_format=fmt, rows=rows)


# ---------------------------------------------------------------------------
# Prokka database converters
# ---------------------------------------------------------------------------


def genbank_to_protein_fasta_prokka(
    genbank_files: Sequence[str],
    output_path: str,
    *,
    id_tag: str | None = None,
    include_hypothetical: bool = False,
    include_pseudo: bool = False,
    min_length: int = 0,
    genetic_code: int | None = None,
    separator: str = "~~~",
    blank: str = "",
    prokka_script_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """GenBank CDS features to a protein FASTA reference (Prokka ``prokka-genbank_to_fasta_db``).

    Each CDS is translated (``/transl_table`` or ``genetic_code``) and written as ``>id EC~~~gene~~~product``. ``id_tag``
    selects the qualifier used as FASTA id (default: first of protein_id, locus_tag, db_xref).
    """
    sources = cli.existing_files(list(genbank_files), "genbank_files")
    target = cli.new_file(output_path, "output_path")
    arguments = ["--sep", separator if separator else "~~~", "--minlen", str(cli.integer(min_length, "min_length", minimum=0))]
    if id_tag is not None:
        arguments += ["--idtag", cli.text(id_tag, "id_tag")]
    if genetic_code is not None:
        arguments += ["--gcode", str(cli.integer(genetic_code, "genetic_code", minimum=1, maximum=33))]
    if blank:
        arguments += ["--blank", cli.text(blank, "blank")]
    arguments.append("--hypo" if cli.flag(include_hypothetical, "include_hypothetical") else "--nohypo")
    arguments.append("--pseudo" if cli.flag(include_pseudo, "include_pseudo") else "--nopseudo")
    arguments += [str(p) for p in sources]
    record = cli.invoke(_PROKKA_DB, arguments, cwd=target.parent, stdout_to=target, explicit=prokka_script_path)
    return cli.finish(record, "genbank_to_protein_fasta_prokka", proteins=sum(1 for line in target.open() if line.startswith(">")))


def translate_abricate_genes_prokka(
    abricate_sequences: str,
    output_path: str,
    *,
    prokka_script_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Translate an ABRicate nucleotide gene FASTA into a Prokka protein reference (``prokka-abricate_to_fasta_db``).

    Input headers must be ABRicate's ``db~~~gene~~~accession~~~description``. The program forces a Met start, removes the
    stop codon and refuses sequences with internal stops or duplicate proteins.
    """
    source = cli.existing_file(abricate_sequences, "abricate_sequences")
    target = cli.new_file(output_path, "output_path")
    record = cli.invoke(_PROKKA_DB, [str(source)], cwd=target.parent, stdout_to=target, explicit=prokka_script_path,
                        executable_name="prokka-abricate_to_fasta_db")
    return cli.finish(record, "translate_abricate_genes_prokka", proteins=sum(1 for line in target.open() if line.startswith(">")))


__all__ = [
    "add_kaiju_taxon_names",
    "annotate_dna_eggnog_mapper",
    "build_bracken_kmer_distribution",
    "check_gff_annotation_amrfinder",
    "check_nucleotide_fasta_amrfinder",
    "combine_bracken_tables",
    "convert_kaiju_to_krona",
    "detect_amr_genes_nucleotide",
    "download_abricate_database",
    "genbank_to_protein_fasta_prokka",
    "map_fastani_conserved_regions",
    "merge_kaiju_outputs",
    "plot_fastani_conserved_regions",
    "reestimate_bracken_report",
    "run_genomad",
    "screen_contigs_abricate",
    "search_dna_translated_diamond",
    "search_kaiju_best_match",
    "summarize_abricate_reports",
    "summarize_kaiju_table",
    "translate_abricate_genes_prokka",
    "type_sequence_type_mlst",
    "type_shiga_toxin_operons",
]

_register_public_verbose(globals(), tuple(__all__))
del os
